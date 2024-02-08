<?php

/**
 * @page analyze_rna
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("analyze_rna", "RNA mapping pipeline using STAR.");
$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);
//optional
$parser->addInfile("system", "Processing system INI file (determined from NGSD via the 'name' by default).", true);
$steps_all = array("ma", "rc", "an", "fu", "db", "plt");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, rc=read counting, an=annotation, fu=fusion detection, db=import into NGSD", true, "ma,rc,an,fu,db,plt");
$parser->addEnum("library_type", "Specify the library type, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addFlag("skip_dedup", "Skip alignment duplication marking (coordinate or UMI based).");
$parser->addFlag("skip_filter_hb", "Do not automatically filter input FASTQ for globin reads for blood samples.");
$parser->addString("out_folder", "Folder where analysis results should be stored. Default is same as in '-folder' (e.g. Sample_xyz/).", true, "default");
$parser->addInt("threads", "The maximum number of threads to use.", true, 5);
$parser->addFlag("skip_dna_reannotation", "Do not automatically start the reannotation of the related DNA sample.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");

extract($parser->parse($argv));

//resolve out_folder
if($out_folder=="default")
{
	$out_folder = $folder;
}
$ngsbits = get_path("ngs-bits");

// create logfile in output folder if no filepath is provided:
if ($parser->getLogFile() == "") $parser->setLogFile($out_folder."/analyze_rna_".date("YmdHis").".log");

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all))
	{
		trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	}
}

//log server name
list($server) = exec2("hostname -f");
$user = exec('whoami');
$parser->log("Executed on server: ".implode(" ", $server)." as ".$user);

//init
$prefix = $out_folder."/".$name;
$sys = load_system($system, $name);
$build = $sys['build'];
$target_file = $sys['target_file'];

//determine genome from build
$genome = get_path("data_folder")."/genomes/STAR/{$build}/";

//determine gtf from build
$gtfFile = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$build);
}

//find FASTQ files
$in_for = glob($folder."/*_R1_001.fastq.gz");
$in_rev = glob($folder."/*_R2_001.fastq.gz");
$in_umi = glob($folder."/*_index_*.fastq.gz");
//set UMI flag
$umi = in_array($sys['umi_type'], ["IDT-UDI-UMI"]) && !empty($in_umi);

if (empty($in_for) && empty($in_rev))
{
	trigger_error("No FASTQ files found!", E_USER_ERROR);
}

$filter_hb = false;
if (db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
	$info = get_processed_sample_info($db, $name, false);
	if (!is_null($info))
	{
		$filter_hb = (strcasecmp($info["tissue"], "blood") == 0) && !$skip_filter_hb;
		if ($filter_hb) trigger_error("Sample type is blood, enabling HB read filter", E_USER_NOTICE);
	}
}

//mapping and QC
$final_bam = $prefix.".bam";
$out_splicing = "{$prefix}_splicing.tsv";
$out_chimeric = "{$prefix}_chimeric.tsv";
$qc_fastq = $prefix."_stats_fastq.qcML";
$qc_map = $prefix."_stats_map.qcML";
if (in_array("ma", $steps))
{
	//check FASTQ quality encoding
	$files = array_merge($in_for, $in_rev);
	foreach ($files as $file)
	{
		list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."FastqFormat", "-in {$file}", true);
		if (!contains($stdout[2], "Sanger"))
		{
			trigger_error("Input file '{$file}' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
		}
	}

	//check that adapters are specified
	if ($sys["adapter1_p5"]=="" || $sys["adapter2_p7"]=="")
	{
		trigger_error("No forward and/or reverse adapter sequence given!\nForward: ".$sys["adapter1_p5"]."\nReverse: ".$sys["adapter2_p7"], E_USER_ERROR);
	}

	//adapter trimming + QC
	$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
	$fastq_trimmed2 = $parser->tempFile("_trimmed.fastq.gz");
	$seqpurge_params = [
		"-out1", $fastq_trimmed1,
		"-out2", $fastq_trimmed2,
		"-a1", $sys["adapter1_p5"],
		"-a2", $sys["adapter2_p7"],
		"-qc", $qc_fastq,
		"-threads", bound($threads, 1, 6),
		"-qcut 0"
	];

	if (in_array($sys['umi_type'], ["IDT-UDI-UMI"]))
	{
		if (empty($in_umi))
		{
			// check if UMIs are present in FastQ header
			list($stdout, $stderr, $return_code) = $parser->exec(get_path("ngs-bits")."FastqCheckUMI", "-in ".$in_for[0]);

			trigger_error("FastqCheckUMI output:\t".implode("\n", $stdout), E_USER_NOTICE);
			if(starts_with($stdout[0], "UMI: true"))
			{
				$umi = true;
				trigger_error("UMIs annotated to the FastQ header will be used.", E_USER_NOTICE);
			}
			else
			{
				trigger_error("No UMI read files found! Aborting analysis.", E_USER_ERROR);
			}

			$seqpurge_params = array_merge($seqpurge_params,
			[
				"-in1", implode(" ", $in_for),
				"-in2", implode(" ", $in_rev)
			]);
		}
		else
		{
			// add barcodes to header
			$merged1_bc = $parser->tempFile("_bc1.fastq.gz");
			$merged2_bc = $parser->tempFile("_bc2.fastq.gz");
			$parser->exec(get_path("ngs-bits")."FastqAddBarcode", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -in_barcode ".implode(" ", $in_umi)." -out1 $merged1_bc -out2 $merged2_bc", true);
			
			$seqpurge_params = array_merge($seqpurge_params,
				[
					"-in1", $merged1_bc,
					"-in2", $merged2_bc
				]);
		}
	}
	else
	{
		$seqpurge_params = array_merge($seqpurge_params,
			[
				"-in1", implode(" ", $in_for),
				"-in2", implode(" ", $in_rev)
			]);
	}



	$parser->exec(get_path("ngs-bits")."SeqPurge", implode(" ", $seqpurge_params), true);

	}
	if (in_array("ma", $steps))
	{
		if ($filter_hb)
		{
			$kraken_tmpdir = $parser->tempFolder("kraken2_filter_hb");
			$filtered = "{$kraken_tmpdir}/filtered#.fastq";
			$filtered1 = "{$kraken_tmpdir}/filtered_1.fastq";
			$filtered2 = "{$kraken_tmpdir}/filtered_2.fastq";
			$kraken_args = [
				"--db", get_path("data_folder")."/dbs/kraken2_filter_hb",
				"--threads", $threads,
				"--output", "-",
				"--paired",
				"--gzip-compressed",
				"--unclassified-out", "{$kraken_tmpdir}/filtered#.fastq",
				$fastq_trimmed1,
				$fastq_trimmed2
			];
			$parser->exec(get_path("kraken2"), implode(" ", $kraken_args));
			$parser->exec("gzip", "-1 {$filtered1}");
			$parser->exec("gzip", "-1 {$filtered2}");

			$fastq_trimmed1 = "{$filtered1}.gz";
			$fastq_trimmed2 = "{$filtered2}.gz";
		}
		//mapping
		$tmp_aligned = $parser->tempFile("aligned.bam");
		$args = array(
			"-out", $tmp_aligned,
			"-threads", $threads,
			"-in1", $fastq_trimmed1,
			"-in2", $fastq_trimmed2,
			"-genome", $genome,
			"-out_splicing", $out_splicing,
			"-out_chimeric", $out_chimeric,
			"--log", $parser->getLogFile()
		);

		//skip duplicate marking if requested, or for UMIs (will be handled later)
		if ($skip_dedup || $umi) $args[] = "-skip_dedup";

		$parser->execTool("NGS/mapping_star.php", implode(" ", $args));

		if ($umi && !$skip_dedup)
		{
			//UMI-based duplicate flagging
			$pipeline = [];
			//UMI-tools group alignments (coordinate sorted)
			$pipeline[] = [get_path("umi_tools"), "group --stdin {$tmp_aligned} --output-bam --log2stderr --paired --no-sort-output --umi-separator=':' --compresslevel 0 --unmapped-reads use"];
			$tmp_for_sorting1 = $parser->tempFile();
			//sort by query name
			$pipeline[] = [get_path("samtools"), "sort -n -T {$tmp_for_sorting1} -m 2G -@ {$threads}"];
			//use fixmate to add mate scores
			$pipeline[] = [get_path("samtools"), "fixmate -@ {$threads} -u -m - -"];
			$tmp_for_sorting2 = $parser->tempFile();
			//sort by coordinate
			$pipeline[] = [get_path("samtools"), "sort -T {$tmp_for_sorting2} -m 2G -@ {$threads}"];
			$tmp_for_markdup = $parser->tempFile();
			//mark duplicates by UG tag (unique group from umi_tools group)
			$pipeline[] = [get_path("samtools"), "markdup -T {$tmp_for_markdup} -c --barcode-tag UG -s -@ {$threads} - {$final_bam}"];
			$parser->execPipeline($pipeline, "umi-tools dedup pipeline");

			//index
			$parser->indexBam($final_bam, $threads);
		}
		else
		{
			$parser->moveFile($tmp_aligned, $final_bam);
			$parser->moveFile("{$tmp_aligned}.bai", "{$final_bam}.bai");
		}

		//mapping QC
		$mappingqc_params = array(
			"-in ".$final_bam,
			"-out ".$qc_map,
			"-ref ".genome_fasta($sys['build']),
			"-build ".ngsbits_build($sys['build'])
		);

		$mappingqc_params[] = (isset($target_file) && $target_file != "") ? "-roi {$target_file}" : "-rna";
		if ($build!="GRCh38") $mappingqc_params[] = "-no_cont";

		$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $mappingqc_params), true);
	}

//read counting
$counts_raw = $prefix."_counts_raw.tsv";
$counts_exon_raw = $prefix."_counts_exon_raw.tsv";
$counts_normalized = $prefix."_counts.tsv";
$counts_exon_normalized = $prefix."_counts_exon.tsv";
$counts_qc = $prefix."_stats_rc.tsv";
$repair_bam = $final_bam;
if (in_array("rc", $steps))
{
	$tmpdir = $parser->tempFolder();
	$repair_bam = "{$tmpdir}/{$name}.bam";
	$parser->exec(dirname(get_path("feature_counts")) . "/utilities/repair",
		"-i {$final_bam} -o {$repair_bam} -T {$threads}", true);

	$args_common = array(
		"-in", $repair_bam,
		"-library_type", $library_type,
		"-gtf_file", $gtfFile,
		"-threads", $threads
	);

	// for UMI libraries, ignore duplicates
	if ($umi) $args_common[] = "-ignore_dup";

	$args = array_merge($args_common, [
		"-out", $counts_raw,
		"-qc_file", $counts_qc
	]);

	// gene-level counting
	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args));

	// exon-level counting
	$args_exon = array_merge($args_common, [
		"-exon_level",
		"-out", $counts_exon_raw
	]);
	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args_exon));

	// read count normalization
	$parser->execTool("NGS/rc_normalize.php", "-in $counts_raw -out $counts_normalized -in_exon $counts_exon_raw -out_exon $counts_exon_normalized");
}

//annotate
$expr = $prefix."_expr.tsv";
$expr_exon = $prefix."_expr_exon.tsv";
$expr_stats = $prefix."_expr.stats.tsv";
$expr_corr = $prefix."_expr.corr.txt";
$junctions = "{$prefix}_splicing.tsv";
$splicing_annot = "{$prefix}_splicing_annot.tsv";
$splicing_bed = "{$prefix}_splicing.bed";
$splicing_gene = "{$prefix}_splicing_gene.tsv";

$reference_tissues = array();
if (in_array("an", $steps))
{
	//annotate gene-level read counts
	$parser->execTool("NGS/rc_annotate.php", "-in $counts_normalized -out $counts_normalized -gtfFile $gtfFile -annotationIds gene_name,gene_biotype");
	//annotate exon-level read counts
	$parser->execTool("NGS/rc_annotate.php", "-in $counts_exon_normalized -out $counts_exon_normalized -gtfFile $gtfFile -annotationIds gene_name,gene_biotype");

	//expression value based on cohort
	if (db_is_enabled("NGSD"))
	{
		$db = DB::getInstance("NGSD");
		$ps_info = get_processed_sample_info($db, $name, false);
		if (!is_null($ps_info))
		{
			if ($build=="GRCh38")
			{
				$cohort_strategy = "RNA_COHORT_GERMLINE";
				$hpa_parameter = "";
				if ($ps_info['is_tumor']) 
				{
					$rna_id = get_processed_sample_id($db, $name);
					$s_id = $db->getValue("SELECT sample_id FROM processed_sample where id=$rna_id");
					$reference_tissues = $db->getValues("SELECT DISTINCT sdi.disease_info FROM sample s LEFT JOIN sample_relations sr ON s.id=sr.sample1_id OR s.id=sr.sample2_id LEFT JOIN sample_disease_info sdi ON sdi.sample_id=sr.sample1_id OR sdi.sample_id=sr.sample2_id WHERE s.id=$s_id AND sdi.type='RNA reference tissue' AND (sr.relation='same sample' OR sr.relation IS NULL)");

					
					$cohort_strategy = "RNA_COHORT_SOMATIC";
					
					if (count($reference_tissues) > 0)
					{
						$hpa_parameter = "-hpa_file ".get_path("data_folder")."/dbs/gene_expression/rna_tissue_consensus_v23.tsv";
					}					
				}
				$parser->exec(get_path("ngs-bits") . "NGSDAnnotateRNA", "-mode genes -update_genes -ps {$name} -cohort_strategy {$cohort_strategy} -in {$counts_normalized} -out {$expr} -corr {$expr_corr} {$hpa_parameter}", true);
				$parser->exec(get_path("ngs-bits") . "NGSDAnnotateRNA", "-mode exons -update_genes -ps {$name} -cohort_strategy {$cohort_strategy} -in {$counts_exon_normalized} -out {$expr_exon}", true);
			}
			
			//remove duplicate exons from file
			$exon_file = file($expr_exon, FILE_IGNORE_NEW_LINES);
			$cache = array();
			$output = array();
			foreach($exon_file as $line)
			{
				if(starts_with($line, "#"))
				{
					$output[] = $line;
				}
				else
				{
					$exon = explode("\t", $line, 1)[0];
					
					if(isset($cache[$exon])) continue;

					$cache[$exon] = true;
					$output[] = $line;
				}
			}
			file_put_contents($expr_exon, implode("\n", $output));
		}
	
		//annotate splice junctions
		if ($build === "GRCh38" && file_exists($junctions))
		{
			$parser->exec("{$ngsbits}SplicingToBed", "-in {$junctions} -report {$splicing_annot} -gene_report {$splicing_gene} -bed {$splicing_bed}", true);
		}
	}
}

// run RNA QC
$qc_rna = "{$prefix}_stats_RNA.qcML";
if (in_array("ma", $steps) || in_array("rc", $steps) || in_array("an", $steps))
{
$args = [
	"-bam", $final_bam,
	"-out", $qc_rna,
	"-ref", genome_fasta($sys['build'])
];
if (ngsbits_build($sys['build']) != "non_human")
{
	$args[] = 	"-housekeeping_genes ".repository_basedir()."/data/gene_lists/housekeeping_genes_".ngsbits_build($sys['build']).".bed";
}
if (isset($target_file) && $target_file != "") 
{
	$args[] = "-roi {$target_file}";
}
if (file_exists($splicing_gene))
{
	$args[] = "-splicing ".$splicing_gene;
}
if (file_exists($expr))
{
	$args[] = "-expression ".$expr;
}
	//run RnaQC
	$parser->exec(get_path("ngs-bits")."RnaQC", implode(" ", $args));
}

//save gene expression plots
$expr_cohort = $prefix."_expr.cohort.tsv";
if (in_array("plt", $steps))
{
	$genelists = glob(repository_basedir()."/data/misc/pathway_genelists/*.txt");
	$all_gene_file = $parser->tempFile(".txt", "genes_");
	$parser->exec("cat", implode(" ", $genelists)." > {$all_gene_file}");
	$cohort_strategy = "RNA_COHORT_GERMLINE";
	$ps_info = get_processed_sample_info($db, $name);
	if ($ps_info['is_tumor']) $cohort_strategy = "RNA_COHORT_SOMATIC";
	$parser->exec(get_path("ngs-bits") . "NGSDExtractRNACohort", "-genes {$all_gene_file} -ps {$name} -sample_expression {$expr} -cohort_strategy {$cohort_strategy} -out {$expr_cohort}", true);

	//create cohort/expr file without comments
	$tmp_cohort = $parser->tempFile(".tsv", "cohort_");
	$parser->exec("egrep", " -v '^##' {$expr_cohort} > {$tmp_cohort}");
	$tmp_expr = $parser->tempFile(".tsv", "expr_");
	$parser->exec("egrep", " -v '^##' {$expr} > {$tmp_expr}");


	$args = [
		"--cohort", $tmp_cohort,
		"--annotation", $tmp_expr,
		"--sample", $name
	];
	$genelists = glob(repository_basedir()."/data/misc/pathway_genelists/*.txt");
	
	//remove old images
	$parser->exec("rm -f", "{$prefix}_expr.*.png");
	
	foreach ($genelists as $genelist)
	{
		$shortname = basename($genelist, ".txt");
		$args_extra = [
			"--genelist", $genelist,
			"--title", "\"Gene Expression: {$shortname}\"",
			"--plot", "{$prefix}_expr.{$shortname}.png",
			"--reference"
		];
		$parser->exec(get_path("python3")." ".repository_basedir()."/src/NGS/rc_plot_expr.py", implode(" ", array_merge($args, $args_extra)), true);
	}
}

//detect fusions
$fusions_arriba_tsv = "{$prefix}_fusions_arriba.tsv";
$fusions_arriba_vcf = "{$prefix}_fusions_arriba.vcf";
$fusions_arriba_bam = "{$prefix}_fusions_arriba.bam";
$fusions_arriba_discarded_tsv = "{$prefix}_fusions_arriba.discarded.tsv";
$fusions_arriba_pdf = "{$prefix}_fusions_arriba.pdf";
$fusions_arriba_pic_dir = "{$prefix}_fusions_arriba_pics";
if (in_array("fu",$steps))
{
	$arriba_args = [
		"-bam", $final_bam,
		"-out_fusions", $fusions_arriba_tsv,
		"-out_vcf", $fusions_arriba_vcf,
		"-out_discarded", $fusions_arriba_discarded_tsv,
		"-out_bam", $fusions_arriba_bam,
		"-out_pdf", $fusions_arriba_pdf,
		"-out_pic_dir", $fusions_arriba_pic_dir,
		"-build", $build,
		"--log", $parser->getLogFile()
	];
	$parser->execTool("NGS/vc_arriba.php", implode(" ", $arriba_args));
}

//import to database
if (in_array("db", $steps))
{
	$parser->execTool("NGS/db_check_gender.php", "-in $final_bam -pid $name --log ".$parser->getLogFile());

	//import QC data
	$qc_file_list = array($qc_fastq, $qc_map);
	if (file_exists($qc_rna)) $qc_file_list[] = $qc_rna;
	$parser->exec("{$ngsbits}/NGSDImportSampleQC", "-ps $name -files ".implode(" ", $qc_file_list)." -force");

	//import expression data
	$args = [
		"-ps", $name,
		"-force"
	];
	if ($build=="GRCh38")
	{
		if (file_exists($expr)) $parser->exec(get_path("ngs-bits")."NGSDImportExpressionData", "-mode genes -expression {$expr} ".implode(" ", $args));
		if (file_exists($expr_exon)) $parser->exec(get_path("ngs-bits")."NGSDImportExpressionData", "-mode exons -expression {$expr_exon} ".implode(" ", $args));
	}	
}


if (! $skip_dna_reannotation && db_is_enabled("NGSD") && (in_array("ma", $steps) || in_array("rc", $steps) || in_array("an", $steps) || in_array("fu", $steps)) && count($reference_tissues) > 0)
{
	$db = DB::getInstance("NGSD", false);
	$ps_info = get_processed_sample_info($db, $name);
	$sample_id = $ps_info["s_id"];

	$existing_relations = $db->executeQuery("SELECT * FROM sample_relations WHERE sample1_id = '{$sample_id}' OR sample2_id = '{$sample_id}'");

	$already_started = array();

	foreach($existing_relations as $relation)
	{
		if ($relation["relation"] != "same sample") continue;
		
		$related_sample_id = $relation["sample1_id"];
		if($relation["sample1_id"] == $sample_id)
		{
			$related_sample_id = $relation["sample2_id"];
		}
		$related_sample_info = $db->executeQuery("SELECT name, tumor FROM sample WHERE id = '$related_sample_id'")[0];
		$related_sample_name = $related_sample_info["name"];
		$is_tumor = $related_sample_info["tumor"] != "0";
		$related_processed_samples = $db->executeQuery("SELECT * FROM processed_sample WHERE sample_id = '$related_sample_id'");

		foreach($related_processed_samples as $rps)
		{
			$rps_name = $related_sample_name."_0".$rps["process_id"];
			$rps_info = get_processed_sample_info($db, $rps_name);

			if (! file_exists($rps_info["ps_folder"]) || in_array($rps, $already_started)) continue;

			$output = array();
			if ($is_tumor)
			{
				//search normal sample:
				$normal_ps_id = $rps["normal_id"];
				$normal_infos = $db->executeQuery("SELECT s.name, ps.process_id FROM processed_sample as ps, sample as s WHERE ps.id = '$normal_ps_id' and ps.sample_id = s.id")[0];
				$normal_ps_name = $normal_infos["name"]."_0".$normal_infos["process_id"];
				$output = $parser->execTool("NGS/db_queue_analysis.php", "-user unknown -type 'somatic' -samples $rps_name $normal_ps_name -info tumor normal -args '-steps an_rna'", false);
			}
			else
			{
				$output = $parser->execTool("NGS/db_queue_analysis.php", "-user unknown -type 'single sample' -samples $rps_name -args '-steps vc -annotation_only'", false);
			}
			
			if (count($output[1])>0)
			{
				trigger_error("Error while queueing reannotation of DNA sample $rps_name: \n".implode("\n",$output[1])."With return value ".$output[2].".", E_USER_WARNING);
			}
			else
			{
				print "Starting reannotation for sample: $rps_name.\n";
			}
			$already_started[] = $rps;
		}
	}
}
?>
