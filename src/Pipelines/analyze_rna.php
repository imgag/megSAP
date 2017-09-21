<?php

/**
	@page analyze_rna
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
$steps_all = array("ma", "rc", "an", "fu", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, rc=read counting, an=annotation, fu=fusion detection, db=import into NGSD", true, implode(",", $steps_all));

$parser->addEnum("library_type", "Specify the library type, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addFlag("no_splicing", "Disable spliced read alignment.");
$parser->addFlag("abra", "Enable indel realignment with ABRA.");

$parser->addString("out_folder", "Folder where analysis results should be stored. Default is same as in '-folder' (e.g. Sample_xyz/).", true, "default");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);

extract($parser->parse($argv));

//resolve out_folder
if($out_folder=="default")
{
	$out_folder = $folder;
}

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

if ((count($in_for) == 0) && (count($in_rev) == 0))
{
	trigger_error("No FASTQ files found!", E_USER_ERROR);
}

$paired = (count($in_rev) != 0);

//mapping and QC
$log_ma  = $prefix."_log1_map.log";
$final_bam = $prefix.".bam";
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

	//adapter trimming + QC (SeqPurge for paired-end, ReadQC+skewer for single-end)
	if ($paired)
	{
		$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
		$fastq_trimmed2 = $parser->tempFile("_trimmed.fastq.gz");
		$seqpurge_params = array(
			"-in1", implode(" ", $in_for),
			"-in2", implode(" ", $in_rev),
			"-out1 {$fastq_trimmed1}",
			"-out2 {$fastq_trimmed2}",
			"-a1", $sys["adapter1_p5"],
			"-a2", $sys["adapter2_p7"],
			"-qc", $qc_fastq,
			"-threads", min($threads, 2), // use at most 2 threads
			"-qcut 0"
			);
		$parser->exec(get_path("ngs-bits")."SeqPurge", implode(" ", $seqpurge_params), true);
	}
	else
	{
		$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 ".implode(" ", $in_for)." -out $qc_fastq", true);

		$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
		
		$skewer_params = array(
			"-x", $sys["adapter1_p5"],
			"-y", $sys["adapter2_p7"],
			"-m any",
			"--threads", min($threads, 2), // use at most 2 threads
			"--stdout",
			"-"
		);
		
		$pipline = array();
		$pipeline[] = array("zcat", implode(" ", $in_for));
		$pipeline[] = array(get_path("skewer"), implode(" ", $skewer_params));
		$pipeline[] = array("gzip", "-1 > {$fastq_trimmed1}");
		$parser->execPipeline($pipeline, "skewer");
	}

	//mapping
	$args = array(
		"-out", $final_bam,
		"-threads", $threads,
		"-in1", $fastq_trimmed1,
		"-genome", $genome,
		"--log", $log_ma
	);

	if ($paired) $args[] = "-in2 $fastq_trimmed2";
	if ($no_splicing) $args[] = "-no_splicing";
	
	if (file_exists($log_ma)) unlink($log_ma);
	$parser->execTool("NGS/mapping_star.php", implode(" ", $args));

	//indel realignment
	if ($abra)
	{
		$junction_file = "{$prefix}_splicing.tsv";
		if (!file_exists($junction_file))
		{
			trigger_error("Could not open junction file '$junction_file' needed for indel realignment. Please re-run mapping step.", E_USER_ERROR);
		}
		
		$abra_out = $parser->tempFile("_abra_realigned.bam");
		
		$abra_params = array(
			"-in", $final_bam,
			"-out", $abra_out,
			"-threads", $threads,
			"-build", $build,
			"-gtf", $gtfFile,
			"-junctions", $junction_file
		);
		if (!$paired) $abra_params[] = "-se";
		if (isset($target_file) && $target_file != "") $abra_params[] = "-roi {$target_file}";
		
		$parser->execTool("NGS/indel_realign_abra.php", implode(" ", $abra_params));

		copy2($abra_out, $final_bam);
		$parser->exec(get_path("samtools")." index", $final_bam, true);
	}

	//mapping QC
	$mappingqc_params = array(
		"-in", $final_bam,
		"-out", $qc_map
	);

	$mappingqc_params[] = (isset($target_file) && $target_file != "") ? "-roi {$target_file}" : "-rna";
	if (!in_array($build, array("GRCh37", "hg19"))) $mappingqc_params[] = "-no_cont";

	$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $mappingqc_params), true);
}

//read counting
$counts_raw = $prefix."_counts_raw.tsv";
$counts_normalized = $prefix."_counts.tsv";
if (in_array("rc", $steps))
{
	$args = array(
		"-in", $final_bam,
		"-out", $counts_raw,
		"-library_type", $library_type,
		"-gtf_file", $gtfFile,
		"-threads", $threads
	);
	
	if (!$paired) $args[] = "-single_end";
	
	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args));

	//normalize read counts
	$parser->execTool("NGS/rc_normalize.php", "-in $counts_raw -out $counts_normalized");
}

//annotate
if (in_array("an", $steps))
{
	$parser->execTool("NGS/rc_annotate.php", "-in $counts_normalized -out $counts_normalized -gtfFile $gtfFile");
}

//detect fusions
if (in_array("fu",$steps))
{
	//path of STAR-Fusion index
	$fusion_index = get_path("data_folder")."/genomes/STAR-Fusion/{$build}";
	
	if (is_dir($fusion_index)) {
	
		//add samtools to path
		putenv("PATH=".dirname(get_path("samtools")).":".getenv("PATH"));

		$fusion_tmp_folder = $parser->tempFolder();
		$chimeric_file = "{$prefix}_chimeric.tsv";
		if (!file_exists($chimeric_file))
		{
			trigger_error("Could not open chimeric file '$chimeric_file' needed for STAR-Fusion. Please re-run mapping step.", E_USER_ERROR);
		}

		//remove header from chimeric file for STAR-Fusion
		$chimeric_file_tmp = $parser->tempFile("_STAR_chimeric.tsv");
		$chimeric_file_f = file($chimeric_file);
		array_shift($chimeric_file_f);
		file_put_contents($chimeric_file_tmp, $chimeric_file_f);

		$starfusion_params = array(
			"--genome_lib_dir", $fusion_index,
			"-J", $chimeric_file_tmp,
			"--output_dir", $fusion_tmp_folder
		);

		$parser->exec(get_path("STAR-Fusion"), implode(" ", $starfusion_params), true);
		$parser->exec("cp", "{$fusion_tmp_folder}/star-fusion.fusion_candidates.final.abridged {$prefix}_var_fusions.tsv", true);
	
	}
	else
	{
		trigger_error("STAR-Fusion index not present. Skipping fusion detection step.", E_USER_WARNING);
	}
}

//import to database
$log_db  = $prefix."_log_db.log";
if (in_array("db", $steps))
{
	if(file_exists($log_db)) unlink($log_db);
	$parser->execTool("NGS/db_check_gender.php", "-in $final_bam -pid $name --log $log_db");

	//import QC data
	$parser->execTool("NGS/db_import_qc.php", "-id $name -files $qc_fastq $qc_map -force --log $log_db");

	//update last analysis date
	updateLastAnalysisDate($name, $final_bam);
}