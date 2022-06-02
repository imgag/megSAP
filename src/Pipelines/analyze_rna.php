<?php

/**
 * @page analyze_rna
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("analyze_rna", "RNA mapping pipeline using STAR.");

//mandatory
$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);

//optional
$parser->addInfile("system", "Processing system INI file (determined from NGSD via the 'name' by default).", true);
$steps_all = array("ma", "rc", "an", "fu", "db", "plt");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, rc=read counting, an=annotation, fu=fusion detection, db=import into NGSD", true, "ma,rc,an,fu,db,plt");

$parser->addEnum("library_type", "Specify the library type, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addFlag("skip_dedup", "Skip alignment duplication marking.");
$parser->addFlag("skip_filter_hb", "Do not automatically filter input FASTQ for globin reads for blood samples.");

$parser->addString("out_folder", "Folder where analysis results should be stored. Default is same as in '-folder' (e.g. Sample_xyz/).", true, "default");
$parser->addInt("threads", "The maximum number of threads to use.", true, 5);

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
		$filter_hb = ($info["tissue"] == "Blood") && !$skip_filter_hb;
		trigger_error("Sample type is blood, enabling HB read filter", E_USER_NOTICE);
	}
}

//mapping and QC
$final_bam = $prefix.".bam";
$before_dedup_bam = $prefix."_before_dedup.bam";
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
			trigger_error("No UMI read files found! Processing fastqs without UMIs.", E_USER_WARNING);
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
		// $parser->exec("pigz", "-p {$threads} {$filtered1}");
		// $parser->exec("pigz", "-p {$threads} {$filtered2}");

		$fastq_trimmed1 = "{$filtered1}.gz";
		$fastq_trimmed2 = "{$filtered2}.gz";
	}
	//mapping
	$args = array(
		"-out", $umi ? $before_dedup_bam : $final_bam,
		"-threads", $threads,
		"-in1", $fastq_trimmed1,
		"-in2", $fastq_trimmed2,
		"-genome", $genome,
		"--log", $parser->getLogFile()
	);

	if ($skip_dedup) $args[] = "-skip_dedup";
	
	$parser->execTool("NGS/mapping_star.php", implode(" ", $args));

	if ($umi)
	{
		//generate $final_bam from $before_dedup_bam

		//barcode correction
		$pipeline = [];
		//UMI-tools dedup
		$pipeline[] = [get_path("umi_tools"), "dedup --stdin {$before_dedup_bam} --out-sam --log2stderr --paired --mapping-quality=3 --no-sort-output --umi-separator=':' --output-stats={$prefix}_umistats"];
		//remove DUP flags
		// TODO replace with samtools view --remove-flags, new samtools version required
		$pipeline[] = [get_path("samtools"), "view --remove-flags DUP -u -b"];
		//sort
		$tmp_for_sorting = $parser->tempFile();
		$pipeline[] = [get_path("samtools"), "sort -T {$tmp_for_sorting} -m 1G -@ ".min($threads, 4)." -o {$final_bam} -"];
		$parser->execPipeline($pipeline, "umi-tools dedup pipeline");

		//index
		$parser->indexBam($final_bam, $threads);
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

	$args = array_merge($args_common, [
		"-out", $counts_raw,
		"-qc_file", $counts_qc
	]);

	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args));

	// exon-level counting
	$args_exon = array_merge($args_common, [
		"-exon_level",
		"-out", $counts_exon_raw
	]);
	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args_exon));

	// read count normalization
	$parser->execTool("NGS/rc_normalize.php", "-in $counts_raw -out $counts_normalized -in_exon $counts_exon_raw -out_exon $counts_exon_normalized");

	// re-run read counting without duplicate alignments
	$counts_nodup = $prefix."_counts_nodup_raw.tsv";
	$args_dup = array_merge($args_common, [
		"-ignore_dup",
		"-out", $counts_nodup
	]);
	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args_dup));
}

//annotate
$expr = $prefix."_expr.tsv";
$expr_cohort = $prefix."_expr.cohort.tsv";
$expr_stats = $prefix."_expr.stats.tsv";
$expr_corr = $prefix."_expr.corr.txt";
$junctions = $umi ? "{$prefix}_before_dedup_splicing.tsv" : "{$prefix}_splicing.tsv";
$splicing_annot = "{$prefix}_splicing_annot.tsv";
$splicing_bed = "{$prefix}_splicing.bed";
$splicing_gene = "{$prefix}_splicing_gene.tsv";
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
			$args = [
				"-name {$name}",
				"-in {$counts_normalized}",
				"-out {$expr}",
				"-cohort {$expr_cohort}",
				"-stats {$expr_stats}",
				"-corr {$expr_corr}"
			];
			//somatic: enable somatic mode, and use RNA reference tissue
			if ($ps_info['is_tumor'])
			{
				$args[] = "-somatic";
				
				
				//somatic case: find HPA reference tissue
				list($s_name) = explode("_", $name);
				$sql = <<<SQL
					SELECT sdi.disease_info FROM sample s
					LEFT JOIN sample_relations sr ON s.id=sr.sample1_id OR s.id=sr.sample2_id
					LEFT JOIN sample_disease_info sdi ON sdi.sample_id=sr.sample1_id OR sdi.sample_id=sr.sample2_id
					WHERE
						s.name='{$s_name}' AND
						sdi.type='RNA reference tissue' AND
						(sr.relation="same sample" OR sr.relation IS NULL)
SQL;
		
				$res = array_unique($db->getValues($sql));
				if (count($res) == 1)
				{
					$rna_ref_tissue = $res[0];
					$args[] = "-hpa_tissue '{$rna_ref_tissue}'";
				}
			}
	
			$parser->execTool("NGS/rc_annotate_expr.php", implode(" ", $args));
		}

		//annotate splice junctions
		if ($build === "GRCh38" && db_is_enabled("NGSD") && file_exists($junctions))
		{
			$parser->exec("{$ngsbits}SplicingToBed", "-in {$junctions} -report {$splicing_annot} -gene_report {$splicing_gene} -bed {$splicing_bed}", true);
		}
	}
}

//save gene expression plots
if (in_array("plt", $steps))
{
	$args = [
		"--cohort", $expr_cohort,
		"--annotation", $expr,
		"--sample", $name
	];
	$genelists = glob(get_path("data_folder") . "/pathway_genelists/*.txt");
	foreach ($genelists as $genelist)
	{
		$shortname = basename($genelist, ".txt");
		$args_extra = [
			"--genelist", $genelist,
			"--title", "\"Gene Expression: {$shortname}\"",
			"--plot", "{$prefix}_expr.{$shortname}.png",
			"--reference"
		];
		$parser->exec("python3 ".repository_basedir()."/src/NGS/rc_plot_expr.py", implode(" ", array_merge($args, $args_extra)), true);
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
		"-bam", $umi ? $before_dedup_bam : $final_bam,
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

// run RNA QC
$qc_rna = "{$prefix}_stats_RNA.qcML";
if (in_array("ma", $steps) || in_array("rc", $steps) || in_array("an", $steps))
{
	$args = [
		"-bam", $umi ? $before_dedup_bam : $final_bam,
		"-housekeeping_genes", repository_basedir()."/data/gene_lists/housekeeping_genes_".ngsbits_build($sys['build']).".bed",
		"-out", $qc_rna,
		"-ref", genome_fasta($sys['build'])
	];

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


//import to database
if (in_array("db", $steps))
{
	$parser->execTool("NGS/db_check_gender.php", "-in $final_bam -pid $name --log ".$parser->getLogFile());

	//import QC data
	$qc_file_list = array($qc_fastq, $qc_map);
	if (file_exists($qc_rna)) $qc_file_list[] = $qc_rna;
	$parser->execTool("NGS/db_import_qc.php", "-id $name -files ".implode(" ", $qc_file_list)." -force --log ".$parser->getLogFile());
}

?>