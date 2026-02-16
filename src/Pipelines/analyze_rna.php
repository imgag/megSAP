<?php

/**
 * @page analyze_rna
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("analyze_rna", "RNA mapping pipeline using STAR.");
$parser->addInfile("folder", "Analysis data folder.", false);
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
$parser->addFlag("skip_gender_check", "Skip check of the detected gender against the entry in NGSD.");

extract($parser->parse($argv));


function getCohortStrategy($ps_info)
{
	$cohort_strategy = "RNA_COHORT_GERMLINE";
	if ($ps_info != null && $ps_info['is_tumor']) $cohort_strategy = "RNA_COHORT_SOMATIC"; 
	return $cohort_strategy;
}

function getCohortSamples($db, $ps_name, $cohort_strategy)
{
	global $parser;
	
	$ps_sys = processingSystem($db, $ps_name);
	
	$ps_allowed_systems ="";
	$rna_allowed_systems = get_path("rna_allowed_systems", false);
	if(is_array($rna_allowed_systems) && array_key_exists($ps_sys, $rna_allowed_systems)) $ps_allowed_systems = $rna_allowed_systems[$ps_sys];
	
	$args = [];
	$args[] = "-ps $ps_name";
	$args[] = "-only_samples";
	$args[] = "-cohort_strategy $cohort_strategy";
	$args[] = "-allowed_systems $ps_allowed_systems";
	list ($stdout, $stderr) = execApptainer("ngs-bits", "NGSDExtractRNACohort", implode(" ", $args));
	
	return $stdout;
}

function processingSystem($db, $sample)
{
	$ps_id = get_processed_sample_id($db, $sample);
	return $db->getValue("SELECT sys.name_short FROM processed_sample as ps LEFT JOIN processing_system as sys ON ps.processing_system_id = sys.id WHERE ps.id = $ps_id");
}

function processingSystems($db, $samples)
{
	$processing_systems = [];
	
	foreach($samples as $sample)
	{
		$processing_systems[$sample] = processingSystem($db, $sample);
	}
	
	return $processing_systems;
}

/**
	@brief returns all 'reference'/helper samples to have a enough samples of each system to remove batch effects
	@systems the short_names of processing systems appearing in the cohort
	@return all samples where sys1 and sys2 are in the given systems
*/
function getBatchCorrectSamples($systems, $cohort_samples)
{
	//base samples to be able to remove batch effects between sys1 and sys2.
	$all_samples = [];
	$all_samples[] = ["name1" => "RNA2412583A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2412583A1_03", "sys2" => "Twist_RNA_Exome", "covar" => "lung"];
	$all_samples[] = ["name1" => "RNA2413360A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2413360A1_03", "sys2" => "Twist_RNA_Exome", "covar" => "cerebral_cortex"];
	$all_samples[] = ["name1" => "RNA2500053A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2500053A1_03", "sys2" => "Twist_RNA_Exome", "covar" => "colon"];
	$all_samples[] = ["name1" => "RNA2500212A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2500212A1_03", "sys2" => "Twist_RNA_Exome", "covar" => "cerebral_cortex"];
	$all_samples[] = ["name1" => "RNA2500587A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2500587A1_03", "sys2" => "Twist_RNA_Exome", "covar" => "adipose_tissue"];
	$all_samples[] = ["name1" => "RNA2500692A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2500692A1_03", "sys2" => "Twist_RNA_Exome", "covar" => "tongue"];
	$all_samples[] = ["name1" => "RNA2501696A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2501696A1_03", "sys2" => "Twist_RNA_Exome", "covar" => "rectum"];
	$all_samples[] = ["name1" => "RNA2501794A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2501794A1_03", "sys2" => "Twist_RNA_Exome", "covar" => "skin"];
	$all_samples[] = ["name1" => "RNA2412583A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2412583A1_04", "sys2" => "Twist_RNA_Exome", "covar" => "lung"];
	$all_samples[] = ["name1" => "RNA2413360A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2413360A1_04", "sys2" => "Twist_RNA_Exome", "covar" => "cerebral_cortex"];
	$all_samples[] = ["name1" => "RNA2500053A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2500053A1_04", "sys2" => "Twist_RNA_Exome", "covar" => "colon"];
	$all_samples[] = ["name1" => "RNA2500212A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2500212A1_04", "sys2" => "Twist_RNA_Exome", "covar" => "cerebral_cortex"];
	$all_samples[] = ["name1" => "RNA2500587A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2500587A1_04", "sys2" => "Twist_RNA_Exome", "covar" => "adipose_tissue"];
	$all_samples[] = ["name1" => "RNA2500692A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2500692A1_04", "sys2" => "Twist_RNA_Exome", "covar" => "tongue"];
	$all_samples[] = ["name1" => "RNA2501696A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2501696A1_04", "sys2" => "Twist_RNA_Exome", "covar" => "rectum"];
	$all_samples[] = ["name1" => "RNA2501794A1_01", "sys1" => "nebRNAU2_qiaRRNA_umi", "name2" => "RNA2501794A1_04", "sys2" => "Twist_RNA_Exome", "covar" => "skin"];

	$samples = [];
	foreach($all_samples as $sample)
	{
		if (in_array($sample["sys1"], $systems) && in_array($sample["sys2"], $systems))
		{
			if (!in_array($sample["name1"], $cohort_samples))
			{
				$samples[] = ["name" => $sample["name1"], "sys" => $sample["sys1"], "covar" => $sample["covar"]];
			}
			if (!in_array($sample["name2"], $cohort_samples))
			{
				$samples[] = ["name" => $sample["name2"], "sys" => $sample["sys2"], "covar" => $sample["covar"]];
			}
		}
	}
	
	return $samples;
}

function getCountFile($ps_name, $exons)
{
	$args = [];
	$args[] = "-ps {$ps_name}";
	$args[] = "-type SAMPLE_FOLDER";
	list ($stdout, $stderr) = execApptainer("ngs-bits", "SamplePath", implode(" ", $args));
	$ps_folder = trim(implode("", $stdout));
		
	if ($exons) return $ps_folder.$ps_name."_counts_exon.tsv";
	return $ps_folder.$ps_name."_counts.tsv";
}

function writeCohortCounts($full_cohort, $ps_name, $out, $exons)
{
	global $parser;
	
	$removed_samples = [];
	$cohort_counts = new Matrix();

	//set names based on diag sample:
	$count_file = getCountFile($ps_name, $exons);
	$parser->log("writeCohortCounts(): given sample $ps_name - countfile: $count_file");
	if (! is_file($count_file)) trigger_error("The count file for $ps_name doesn't exist: $count_file!", E_USER_ERROR);
	
	
	$sample_counts = Matrix::fromTSV($count_file);
	$ident_idx = $sample_counts->getColumnIndex("gene_id");
	if ($exons) $ident_idx = $sample_counts->getColumnIndex("gene_exon");
	
	$names = $sample_counts->getCol($ident_idx);
	$cohort_counts->addCol($names, "Identifier");
	
	//gather counts
	foreach($full_cohort as $sample)
	{
		$count_file = getCountFile($sample, $exons);
		
		if (! is_file($count_file))
		{	
			$removed_samples[] = $sample;
			trigger_error("Skipped counts of $sample because file doesn't exist: $count_file\n", E_USER_WARNING);
			continue;
		}
		
		$sample_counts = Matrix::fromTSV($count_file);
		
		$raw_idx = $sample_counts->getColumnIndex("raw");
		$raw_counts = $sample_counts->getCol($raw_idx);
		
		//check count and correct order!
		if (count($raw_counts) == $cohort_counts->rows() && count(array_diff_assoc($names, $sample_counts->getCol($ident_idx))) == 0)
		{
			$cohort_counts->addCol($raw_counts, $sample);
		}
		else
		{
			$removed_samples[] = $sample;
			if (count($raw_counts) != $cohort_counts->rows())
			{
				trigger_error("Skipped counts of $sample because of differing row count: $count_file\n", E_USER_WARNING);
			} 
			else
			{
				trigger_error("Skipped counts of $sample because of differing order of counts: $count_file\n", E_USER_WARNING);
			}
			continue;
		}
	}
	
	$tmp_count_file = $parser->tempFile("_rna_cohort_counts.tsv");
	$cohort_counts->toTSV($tmp_count_file);
	
	$parser->execApptainer("htslib", "bgzip", "-c $tmp_count_file > $out", [], [dirname($out)]);
	
	return $removed_samples;
}

function generateBatchCorrectData($cohort_samples, $ps_name, $processing_systems, $raw_cohort_counts_file, $exons, $batch_correct_samples=[])
{
	$full_cohort = [];
	$full_batches = [];
	
	//determine the covariant of the cohort (tissue) in comparision to added samples.
	$cohort_covar = "other";
	foreach($batch_correct_samples as $ref_sample)
	{
		if (in_array($ref_sample["name"], $cohort_samples))
		{
			$cohort_covar = $ref_sample["covar"];
		}
	}
	
	$full_cohort = $cohort_samples;
	$full_batches = $processing_systems;
	$full_covar = array_fill_keys($cohort_samples, $cohort_covar);
	foreach($batch_correct_samples as $ref_sample)
	{
		$sample = $ref_sample["name"];
		if (! in_array($sample, $full_cohort)) $full_cohort[] = $sample;
		$full_covar[$sample] = $ref_sample["covar"];
		$full_batches[$sample] = $ref_sample["sys"];
	}
	
	//write raw cohort count files
	$removed_samples = writeCohortCounts($full_cohort, $ps_name, $raw_cohort_counts_file, $exons);
	
	//remove samples that couldn't be used in the counts table (wrong number/order of counts -> not comparable because of differing analysis)
	$full_cohort = array_diff($full_cohort, $removed_samples);
	foreach($removed_samples as $remove)
	{
		unset($full_batches[$remove]);
		unset($full_covar[$remove]);
	}
	
	return [$full_cohort, $full_batches, $full_covar];
}

function remove_helper_samples($cohort_count_file, $batch_correct_samples)
{
	global $parser;
	if (! count($batch_correct_samples) > 0) return;
	
	$tmp_cohort_count_file = $parser->tempFile("_analyze_rna_cohort_count_rm_helper.tsv");
	
	$cohort_mat = Matrix::fromTSV($cohort_count_file);
	foreach($batch_correct_samples as $helper_sample)
	{
		$cohort_mat->removeColByName($helper_sample["name"]);
	}
	
	$cohort_mat->toTSV($tmp_cohort_count_file);
	//overwrite old cohort_count file with cleaned version
	$parser->execApptainer("htslib", "bgzip", "-c $tmp_cohort_count_file > $cohort_count_file", [], [dirname($cohort_count_file)]);
}

function writeTSV($array, $out_file, $headers)
{
	$lines = [];
	$lines[] = implode("\t", $headers);
	foreach($array as $sample=>$value)
	{
		$lines[] = $sample."\t".$value;
	}
	
	file_put_contents($out_file, implode("\n", $lines));
}


#//////////////////////////////////
#-----------  MAIN  ---------------
#//////////////////////////////////


//resolve out_folder
if($out_folder=="default")
{
	$out_folder = $folder;
}

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($out_folder."/analyze_rna_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all))
	{
		trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	}
}

//init
$prefix = $out_folder."/".$name;
$sys = load_system($system, $name);
$build = $sys['build'];
$target_file = $sys['target_file'];

//determine genome from build
$genome = get_path("data_folder")."/genomes/STAR/{$build}/";

//determine gtf from build
$gtf_file = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

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
		list($stdout, $stderr) = $parser->execApptainer("ngs-bits", "FastqFormat", "-in {$file}", [$file]);
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
			list($stdout, $stderr, $return_code) = $parser->execApptainer("ngs-bits", "FastqCheckUMI", "-in ".$in_for[0], [$in_for[0]]);

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
			$parser->execApptainer("ngs-bits", "FastqAddBarcode", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -in_barcode ".implode(" ", $in_umi)." -out1 $merged1_bc -out2 $merged2_bc", [$folder]);
			
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


	$parser->execApptainer("ngs-bits", "SeqPurge", implode(" ", $seqpurge_params), [$folder], [$out_folder]);

	//cut off Twist UMIs
	if (in_array($sys['umi_type'], ["Twist"]))
	{
		$cut1 = 5;
		$cut2 = 5;
		$end_trim = 7;
		
		$fastq_trimmed1_bc = $parser->tempFile("_trimmed1_bc.fastq.gz");
		$fastq_trimmed2_bc = $parser->tempFile("_trimmed2_bc.fastq.gz");
		$parser->execApptainer("ngs-bits", "FastqExtractUMI", "-in1 {$fastq_trimmed1} -in2 {$fastq_trimmed2} -out1 {$fastq_trimmed1_bc} -out2 {$fastq_trimmed2_bc} -cut1 {$cut1} -cut2 {$cut2}");

		$parser->deleteTempFile($fastq_trimmed1);
		$parser->deleteTempFile($fastq_trimmed2);

		//special handling of Twist UMIs: cut off first 2bp and last 7bp.
		//trim all ends by 7 bases to remove UMIs from reads with an insert size < read_length
		$fastq_trimmed1_bc_cut = $parser->tempFile("_trimmed1_bc_cut.fastq.gz");
		$fastq_trimmed2_bc_cut = $parser->tempFile("_trimmed2_bc_cut.fastq.gz");
		$parser->execApptainer("ngs-bits", "FastqTrim", "-in {$fastq_trimmed1_bc} -out {$fastq_trimmed1_bc_cut} -start 2 -end $end_trim");
		$parser->execApptainer("ngs-bits", "FastqTrim", "-in {$fastq_trimmed2_bc} -out {$fastq_trimmed2_bc_cut} -start 2 -end $end_trim");
		
		$parser->deleteTempFile($fastq_trimmed1_bc);
		$parser->deleteTempFile($fastq_trimmed2_bc);
		
		$fastq_trimmed1 = $fastq_trimmed1_bc_cut;
		$fastq_trimmed2 = $fastq_trimmed2_bc_cut;
	}


	if ($filter_hb)
	{
		$kraken_tmpdir = $parser->tempFolder("kraken2_filter_hb");
		$filtered = "{$kraken_tmpdir}/filtered#.fastq";
		$filtered1 = "{$kraken_tmpdir}/filtered_1.fastq";
		$filtered2 = "{$kraken_tmpdir}/filtered_2.fastq";

		$kraken_args = [
			"-threads {$threads}",
			"-output -",
			"-paired",
			"-gzip_compressed",
			"-unclassified_out {$kraken_tmpdir}/filtered#.fastq",
			"-in $fastq_trimmed1 $fastq_trimmed2"
		];

		$parser->execTool("Tools/filter_kraken2.php", implode(" ", $kraken_args));

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

	$parser->execTool("Tools/mapping_star.php", implode(" ", $args));

	if ($umi && !$skip_dedup)
	{
		//UMI-based duplicate flagging
		$pipeline = [];
		//UMI-tools group alignments (coordinate sorted)
		$pipeline[] = ["", $parser->execApptainer("umi-tools", "umi_tools group", "--stdin {$tmp_aligned} --output-bam --log2stderr --paired --no-sort-output --umi-separator=':' --compresslevel 0 --unmapped-reads use", [], [], true)];
		$tmp_for_sorting1 = $parser->tempFile();
		//sort by query name
		$pipeline[] = ["", $parser->execApptainer("samtools", "samtools sort", "-n -T {$tmp_for_sorting1} -m 2G -@ {$threads}", [], [], true)];
		//use fixmate to add mate scores
		$pipeline[] = ["", $parser->execApptainer("samtools", "samtools fixmate", "-@ {$threads} -u -m - -", [], [], true)];
		$tmp_for_sorting2 = $parser->tempFile();
		//sort by coordinate
		$pipeline[] = ["", $parser->execApptainer("samtools", "samtools sort", "-T {$tmp_for_sorting2} -m 2G -@ {$threads}", [], [], true)];
		$tmp_for_markdup = $parser->tempFile();
		//mark duplicates by UG tag (unique group from umi_tools group)
		$pipeline[] = ["", $parser->execApptainer("samtools", "samtools markdup", "-T {$tmp_for_markdup} -c --barcode-tag UG -s -@ {$threads} - {$final_bam}", [], [$out_folder], true)];
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
	$in_files = array();
	$mappingqc_params = array(
		"-in ".$final_bam,
		"-out ".$qc_map,
		"-ref ".genome_fasta($sys['build']),
		"-build ".ngsbits_build($sys['build'])
	);
	$in_files[] = $out_folder;
	$in_files[] = genome_fasta($sys['build']);

	if (isset($target_file) && $target_file != "")
	{
		$mappingqc_params[] = "-roi ".realpath($target_file);
		$in_files[] = realpath($target_file);

	}
	else
	{
		$mappingqc_params[] = "-rna";
	}

	if ($build!="GRCh38") $mappingqc_params[] = "-no_cont";

	$parser->execApptainer("ngs-bits", "MappingQC", implode(" ", $mappingqc_params), $in_files);
}

//read counting
$counts_raw = $prefix."_counts_raw.tsv";
$counts_exon_raw = $prefix."_counts_exon_raw.tsv";
$counts_normalized = $prefix."_counts.tsv";
$counts_exon_normalized = $prefix."_counts_exon.tsv";
$counts_qc = $prefix."_stats_rc.tsv";
$repair_bam = $final_bam;
//batch correction files:
$cohort_raw_counts_genes_file = $prefix."_cohort_counts_raw.tsv.gz";
$cohort_raw_counts_exons_file = $prefix."_cohort_counts_exons_raw.tsv.gz";
$cohort_counts_genes_file = $prefix."_cohort_counts.tsv.gz";
$cohort_counts_exons_file = $prefix."_cohort_counts_exons.tsv.gz";
$cohort_counts_normalized_genes_file = $prefix."_cohort_counts_norm.tsv.gz";
$cohort_counts_normalized_exons_file = $prefix."_cohort_counts_exons_norm.tsv.gz";

$uncorrected_counts_genes_normalized = $prefix."_counts_uncorrected.tsv";
$uncorrected_counts_exons_normalized = $prefix."_counts_exon_uncorrected.tsv";
if (in_array("rc", $steps))
{
	$tmpdir = $parser->tempFolder();
	$repair_bam = "{$tmpdir}/{$name}.bam";
	$parser->execApptainer("subread", "repair", "-i {$final_bam} -o {$repair_bam} -T {$threads}", [], [$out_folder]);


	$args_common = array(
		"-in", $repair_bam,
		"-library_type", $library_type,
		"-gtf_file", $gtf_file,
		"-threads", $threads
	);

	// for UMI libraries, ignore duplicates
	if ($umi) $args_common[] = "-ignore_dup";

	$args = array_merge($args_common, [
		"-out", $counts_raw,
		"-qc_file", $counts_qc
	]);

	// gene-level counting
	$parser->execTool("Tools/rc_featurecounts.php", implode(" ", $args));

	// exon-level counting
	$args_exon = array_merge($args_common, [
		"-exon_level",
		"-out", $counts_exon_raw
	]);
	$parser->execTool("Tools/rc_featurecounts.php", implode(" ", $args_exon));
	
	// read count normalization
	$parser->execTool("Tools/rc_normalize.php", "-in $counts_raw -out $counts_normalized -in_exon $counts_exon_raw -out_exon $counts_exon_normalized");
	
	//remove uncorrected files from previous analysis to not interfere with current analysis
	if(is_file($uncorrected_counts_genes_normalized)) exec("rm $uncorrected_counts_genes_normalized");
	if(is_file($uncorrected_counts_exons_normalized)) exec("rm $uncorrected_counts_exons_normalized");
	
	//check if the settings for this processing system allows multiple processing systems in the cohort
	$sys_short_name = $sys["name_short"];
	$ps_allowed_systems ="";
	$rna_allowed_systems = get_path("rna_allowed_systems", false);
	if(is_array($rna_allowed_systems) && array_key_exists($sys_short_name, $rna_allowed_systems)) $ps_allowed_systems = $rna_allowed_systems[$sys_short_name];
	
	if (db_is_enabled("NGSD") && $ps_allowed_systems != "")
	{
		$db = DB::getInstance("NGSD");
		$ps_info = get_processed_sample_info($db, $name, false);
		
		////check if cohort actally differing processing systems
		$cohort_strategy = getCohortStrategy($ps_info);
		$cohort_samples = getCohortSamples($db, $name, $cohort_strategy);
		
		if (! in_array($name, $cohort_samples))
		{
			array_unshift($cohort_samples, $name);
		}
		
		$processing_systems = processingSystems($db, $cohort_samples);
		$processing_system_counts = array_count_values($processing_systems);
		
		//more than one processing system -> remove batch effects
		if (count($processing_system_counts) > 1)
		{
			$add_ref_samples = False;
			foreach(array_keys($processing_system_counts) as $system)
			{
				if ($processing_system_counts[$system] < 10)
				{
					$add_ref_samples = True;
					break;
				}
			}
			
			foreach (["gene", "exon"] as $type)
			{
				$raw_cohort_counts_file = $cohort_raw_counts_genes_file;
				$corrected_cohort_counts_file = $cohort_counts_genes_file;
				$cohort_counts_normalized = $cohort_counts_normalized_genes_file;
				$base_counts_normalized = $counts_normalized;
				$uncorrected = $uncorrected_counts_genes_normalized;
				$raw_count_file = $counts_raw;
				$exons = False;
				
				if ($type == "exon")
				{
					$raw_cohort_counts_file = $cohort_raw_counts_exons_file;
					$corrected_cohort_counts_file = $cohort_counts_exons_file;
					$cohort_counts_normalized = $cohort_counts_normalized_exons_file;
					$base_counts_normalized = $counts_exon_normalized;
					$uncorrected = $uncorrected_counts_exons_normalized;
					$raw_count_file = $counts_exon_raw;
					$exons = True;
				}
				
				//generate necessary data
				$batch_correct_samples = [];
				if ($add_ref_samples)
				{
					//List of samples with keys: [[name, sys, covar], ...]
					$batch_correct_samples = getBatchCorrectSamples(array_keys($processing_system_counts), $cohort_samples);
					list($full_cohort, $full_batches, $full_covar) = generateBatchCorrectData($cohort_samples, $name, $processing_systems, $raw_cohort_counts_file, $exons, $batch_correct_samples);
				}
				else
				{
					list($full_cohort, $full_batches, $full_covar) = generateBatchCorrectData($cohort_samples, $name, $processing_systems, $raw_cohort_counts_file, $exons);
				}
				
				//write tmp files:
				$tmp_batch_file = $parser->tempFile("_batch.tsv");
				$tmp_covar_file = $parser->tempFile("_covar.tsv");
				
				writeTSV($full_batches, $tmp_batch_file, ["SAMPLE", "BATCH"]);
				writeTSV($full_covar, $tmp_covar_file, ["SAMPLE", "COVARIANT"]);
				
				//rename uncorrected counts:
				exec("mv $base_counts_normalized $uncorrected");
				
				$files = [$tmp_batch_file, $tmp_covar_file, $raw_cohort_counts_file, $corrected_cohort_counts_file, repository_basedir()."/src/Tools/rc_batch_correct.py"];
				
				$args = [];
				$args[] = "-i $raw_cohort_counts_file";
				$args[] = "-b $tmp_batch_file";
				if (! $exons) $args[] = "-c $tmp_covar_file";
				$args[] = "-o $corrected_cohort_counts_file";
				$parser->execApptainer("python", "python3", repository_basedir()."/src/Tools/rc_batch_correct.py ".implode(" ", $args), $files);
				
				//remove columns of samples added only to improve batch correction:
				remove_helper_samples($corrected_cohort_counts_file, $batch_correct_samples);
				
				$args = [];
				$args[] = "-ps_name $name";
				$args[] = "-in_cohort $corrected_cohort_counts_file";
				$args[] = "-uncorrected_counts $raw_count_file";				
				$args[] = "-out_cohort $cohort_counts_normalized";
				$args[] = "-out_sample $base_counts_normalized";
				if ($exons)
				{
					$args[] = "-normalized_genes $counts_normalized";
					$parser->execTool("Tools/rc_normalize_cohort_exons.php", implode(" ", $args));
				}
				else
				{
					$parser->execTool("Tools/rc_normalize_cohort_genes.php", implode(" ", $args));
				}
			}
		}
	}	
}

//annotate
$expr = $prefix."_expr.tsv";
$expr_exon = $prefix."_expr_exon.tsv";
$expr_corr = $prefix."_expr.corr.txt";
$junctions = "{$prefix}_splicing.tsv";
$splicing_annot = "{$prefix}_splicing_annot.tsv";
$splicing_bed = "{$prefix}_splicing.bed";
$splicing_gene = "{$prefix}_splicing_gene.tsv";

$uncorrected_expr_genes = $prefix."_expr_uncorrected.tsv";
$uncorrected_expr_exons = $prefix."_expr_exon_uncorrected.tsv";
$uncorrected_expr_corr = $prefix."_expr_uncorrected.corr.txt";

$reference_tissues = array();
if (in_array("an", $steps))
{
	//annotate gene-level read counts
	$parser->execTool("Tools/rc_annotate.php", "-in $counts_normalized -out $counts_normalized -gtf_file $gtf_file -annotationIds gene_name,gene_biotype");
	//annotate exon-level read counts
	$parser->execTool("Tools/rc_annotate.php", "-in $counts_exon_normalized -out $counts_exon_normalized -gtf_file $gtf_file -annotationIds gene_name,gene_biotype");
	
	if (is_file($uncorrected_counts_genes_normalized))
	{
		$parser->execTool("Tools/rc_annotate.php", "-in $uncorrected_counts_genes_normalized -out $uncorrected_counts_genes_normalized -gtf_file $gtf_file -annotationIds gene_name,gene_biotype");
	}
	
	if (is_file($uncorrected_counts_exons_normalized))
	{
		$parser->execTool("Tools/rc_annotate.php", "-in $uncorrected_counts_exons_normalized -out $uncorrected_counts_exons_normalized -gtf_file $gtf_file -annotationIds gene_name,gene_biotype");
	}
	
	//expression value based on cohort
	if (db_is_enabled("NGSD"))
	{
		$in_files = array();
		$db = DB::getInstance("NGSD");
		$ps_info = get_processed_sample_info($db, $name, false);
		if (!is_null($ps_info))
		{
			if ($build=="GRCh38")
			{
				$cohort_strategy = getCohortStrategy($ps_info);
				$hpa_parameter = "";
				if ($ps_info['is_tumor']) 
				{
					$rna_id = get_processed_sample_id($db, $name);
					$s_id = $db->getValue("SELECT sample_id FROM processed_sample where id=$rna_id");
					$reference_tissues = $db->getValues("SELECT DISTINCT sdi.disease_info FROM sample s LEFT JOIN sample_relations sr ON s.id=sr.sample1_id OR s.id=sr.sample2_id LEFT JOIN sample_disease_info sdi ON sdi.sample_id=sr.sample1_id OR sdi.sample_id=sr.sample2_id WHERE s.id=$s_id AND sdi.type='RNA reference tissue' AND (sr.relation='same sample' OR sr.relation IS NULL)");
					
					if (count($reference_tissues) > 0)
					{
						$hpa_parameter = "-hpa_file ".get_path("data_folder")."/dbs/gene_expression/rna_tissue_consensus_v23.tsv";
						$in_files[] = get_path("data_folder")."/dbs/gene_expression/rna_tissue_consensus_v23.tsv";
					}					
				}
				$in_files[] = $out_folder;
				
				if (is_file($uncorrected_counts_genes_normalized))
				{
					$in_files[] = $cohort_counts_normalized_genes_file;
					$parser->execApptainer("ngs-bits", "NGSDAnnotateRNA", "-mode genes -cohort_data $cohort_counts_normalized_genes_file -update_genes -ps {$name} -cohort_strategy {$cohort_strategy} -in {$counts_normalized} -out {$expr} -corr {$expr_corr} {$hpa_parameter}", $in_files);
					$parser->execApptainer("ngs-bits", "NGSDAnnotateRNA", "-mode genes -update_genes -ps {$name} -cohort_strategy {$cohort_strategy} -in {$uncorrected_counts_genes_normalized} -out {$uncorrected_expr_genes} -corr {$uncorrected_expr_corr} {$hpa_parameter}", $in_files);
				}
				else
				{
					$parser->execApptainer("ngs-bits", "NGSDAnnotateRNA", "-mode genes -update_genes -ps {$name} -cohort_strategy {$cohort_strategy} -in {$counts_normalized} -out {$expr} -corr {$expr_corr} {$hpa_parameter}", $in_files);
				}
				
				if (is_file($uncorrected_counts_exons_normalized))
				{
					$parser->execApptainer("ngs-bits", "NGSDAnnotateRNA", "-mode exons -cohort_data $cohort_counts_normalized_exons_file -update_genes -ps {$name} -cohort_strategy {$cohort_strategy} -in {$counts_exon_normalized} -out {$expr_exon}", [$out_folder, $cohort_counts_normalized_genes_file]);
					$parser->execApptainer("ngs-bits", "NGSDAnnotateRNA", "-mode exons -update_genes -ps {$name} -cohort_strategy {$cohort_strategy} -in {$uncorrected_counts_exons_normalized} -out {$uncorrected_expr_exons}", [$out_folder]);
				} 
				else
				{
					$parser->execApptainer("ngs-bits", "NGSDAnnotateRNA", "-mode exons -update_genes -ps {$name} -cohort_strategy {$cohort_strategy} -in {$counts_exon_normalized} -out {$expr_exon}", [$out_folder]);
				}
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
			$parser->execApptainer("ngs-bits", "SplicingToBed", "-in {$junctions} -report {$splicing_annot} -gene_report {$splicing_gene} -bed {$splicing_bed}", [$out_folder]);
		}
	}
}

// run RNA QC
$qc_rna = "{$prefix}_stats_RNA.qcML";
if (in_array("ma", $steps) || in_array("rc", $steps) || in_array("an", $steps))
{
	$in_files = array();
	$args = [
		"-bam", $final_bam,
		"-out", $qc_rna,
		"-ref", genome_fasta($sys['build'])
	];
	$in_files[] = $out_folder;
	$in_files[] = genome_fasta($sys['build']);
	if (ngsbits_build($sys['build']) != "non_human")
	{
		$args[] = 	"-housekeeping_genes ".repository_basedir()."/data/gene_lists/housekeeping_genes_".ngsbits_build($sys['build']).".bed";
		$in_files[] = repository_basedir()."/data/gene_lists/housekeeping_genes_".ngsbits_build($sys['build']).".bed";
	}
	if (isset($target_file) && $target_file != "") 
	{
		$args[] = "-roi {$target_file}";
		$in_files[] = $target_file;
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
		$parser->execApptainer("ngs-bits", "RnaQC", implode(" ", $args), $in_files);
}

//save gene expression plots
$expr_cohort = $prefix."_expr.cohort.tsv";
if (in_array("plt", $steps))
{
	$genelists = glob(repository_basedir()."/data/misc/pathway_genelists/*.txt");
	$all_gene_file = $parser->tempFile(".txt", "genes_");
	$parser->exec("cat", implode(" ", $genelists)." > {$all_gene_file}");
	
	$tmp_expr = $parser->tempFile(".tsv", "expr_");
	$parser->exec("egrep", " -v '^##' {$expr} > {$tmp_expr}");
	
	$tmp_cohort = $parser->tempFile(".tsv", "cohort_");
	if (! file_exists($cohort_counts_normalized_genes_file))
	{
		$ps_info = get_processed_sample_info($db, $name);
		$cohort_strategy = getCohortStrategy($ps_info);
		$parser->execApptainer("ngs-bits", "NGSDExtractRNACohort", "-genes {$all_gene_file} -ps {$name} -sample_expression {$expr} -cohort_strategy {$cohort_strategy} -out {$expr_cohort}", [$out_folder]);
		//create cohort/expr file without comments
		$parser->exec("egrep", " -v '^##' {$expr_cohort} > {$tmp_cohort}");
		
	}
	else
	{
		$tmp_for_cut = $parser->tempFile(".tsv", "cohort_before_cut");
		
		
		$gene_data = $db->executeQuery("SELECT symbol,ensembl_id FROM bioinf_ngsd.gene WHERE ensembl_id is not null");
		
		$gene2ensembl = [];
		foreach($gene_data as $row)
		{
			$gene2ensembl[$row["symbol"]] = $row["ensembl_id"];
		}
		
		$ensemble_ids = array();
		foreach(file($all_gene_file) as $gene)
		{
			$gene = trim($gene);
			if (array_key_exists($gene, $gene2ensembl))
			{
				$ensemble_ids[] = $gene2ensembl[$gene];
			}
			else
			{
				trigger_error("No ensembl ID found for gene {$gene}!", E_USER_ERROR);
			}
		}
		
		list($stdout, $stderr, $return) = $parser->exec("zgrep", "'#gene_id' $cohort_counts_normalized_genes_file");
		$header = str_replace("_tpm", "", trim($stdout[0]));
		file_put_contents($tmp_for_cut, $header."\n");
		
		exec("zgrep '".implode("\|", $ensemble_ids)."' $cohort_counts_normalized_genes_file >> $tmp_for_cut");
		#remove length column:
		exec("cut --complement -f 2 $tmp_for_cut > $tmp_cohort");
		
	}


	$args = [
		"--cohort", $tmp_cohort,
		"--annotation", $tmp_expr,
		"--sample", $name
	];
	
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

		$files = [
			repository_basedir()."/src/Tools/rc_plot_expr.py",
			dirname($genelist),
			$out_folder			
		];

		$parser->execApptainer("python", "python3", repository_basedir()."/src/Tools/rc_plot_expr.py ".implode(" ", array_merge($args, $args_extra)), $files);
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
		"-threads", $threads,
		"--log", $parser->getLogFile(),
		
	];
	$parser->execTool("Tools/vc_arriba.php", implode(" ", $arriba_args));
}

//import to database
if (in_array("db", $steps))
{
	if (! $skip_gender_check)
	{
		$parser->execTool("Tools/db_check_gender.php", "-in $final_bam -pid $name --log ".$parser->getLogFile());
	}
	
	//import QC data
	$qc_file_list = array($qc_fastq, $qc_map);
	if (file_exists($qc_rna)) $qc_file_list[] = $qc_rna;
	$parser->execApptainer("ngs-bits", "NGSDImportSampleQC", "-ps $name -files ".implode(" ", $qc_file_list)." -force", [$out_folder]);

	//import expression data
	$args = [
		"-ps", $name,
		"-force"
	];
	if ($build=="GRCh38")
	{
		//if the NOT batch corrected files exist import those to keep the data in NGSD comparable within the processing systems
		if (file_exists($uncorrected_expr_genes))
		{
			$parser->execApptainer("ngs-bits", "NGSDImportExpressionData", "-mode genes -expression {$uncorrected_expr_genes} ".implode(" ", $args), [$out_folder]);
		}
		else if (file_exists($expr))
		{
			$parser->execApptainer("ngs-bits", "NGSDImportExpressionData", "-mode genes -expression {$expr} ".implode(" ", $args), [$out_folder]);
		}
		
		if (file_exists($uncorrected_expr_exons))
		{
			$parser->execApptainer("ngs-bits", "NGSDImportExpressionData", "-mode exons -expression {$uncorrected_expr_exons} ".implode(" ", $args), [$out_folder]);
		}
		else if (file_exists($expr_exon))
		{
			$parser->execApptainer("ngs-bits", "NGSDImportExpressionData", "-mode exons -expression {$expr_exon} ".implode(" ", $args), [$out_folder]);
		}
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
				$output = $parser->execTool("Tools/db_queue_analysis.php", "-user unknown -type 'somatic' -samples $rps_name $normal_ps_name -info tumor normal -args '-steps an_rna'", false);
			}
			else
			{
				$output = $parser->execTool("Tools/db_queue_analysis.php", "-user unknown -type 'single sample' -samples $rps_name -args '-steps vc -annotation_only'", false);
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
