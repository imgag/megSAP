<?php

/**
	@page queue_rna_sample
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";


error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("queue_rna_sample", "Queue RNA sample for analysis.");

// mandatory arguments
$parser->addString("sample",  "Processed sample identifier.", false);

// optional arguments
$steps_all = array("ma", "rc", "an", "fu", "db");
$parser->addString("steps", "Comma-separated list of processing steps to perform.", true, "ma,rc,an,fu,db");
$parser->addInfile("system",  "Processing system file (determined from database by default).", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 8);
$parser->addString("out_folder", "Output folder for analysis results, defaults to processed sample folder.", true, "");
$parser->addString("add_params", "Additional parameters for the anaylze_rna pipeline.", true, "");
$parser->addFlag("dryrun", "Perform dry-run and only print the command to be queued.");
$parser->addFlag("noqueue", "Do not submit jobs to SGE.");
extract($parser->parse($argv));

// NGSD information
$info = get_processed_sample_info($sample);
$folder = $info['ps_folder'];
$processed_sample = $info['ps_name'];

if(!is_dir($folder))
{
	trigger_error("Could not find sample folder '$sample_folder'!", E_USER_ERROR);
}

if ($out_folder == "") {
	$out_folder = $folder;
}

if (isset($system)) {
	$system_file = $system;
}

// processing system
$sys = load_system($system, $processed_sample);

// input FASTQ files
$in_for = glob($folder."/*_R1_001.fastq.gz");
$in_rev = glob($folder."/*_R2_001.fastq.gz");

if ((count($in_for) == 0) && (count($in_rev) == 0)) {
	trigger_error("No FASTQ files found!", E_USER_ERROR);
}

$paired = (count($in_rev) != 0);

$in_for_s = implode(" ", $in_for);
$in_rev_s = implode(" ", $in_rev);

// start analyze_rna
$analyze_params = array();
$analyze_params[] = "-in_for $in_for_s ";
if ($paired) {
	$analyze_params[] = "-in_rev $in_rev_s";
}
if (isset($system_file)) {
	$analyze_params[] = "-system $system_file";
}
$analyze_params[] = "-out_folder $out_folder -out_name $processed_sample";
$analyze_params[] = "-steps $steps -threads $threads -dedup";
if (isset($sys['stranded']) && $sys['stranded']==1) {
	$analyze_params[] = "-stranded";
}
if (isset($sys['gtfAttribute'])) {
	$analyze_params[] = "-gtfAttribute ".$sys['gtfAttribute'];
}
if (isset($sys['featureType'])) {
	$analyze_params[] = "-featureType ".$sys['featureType'];
}
$analyze_params[] = "-downstream chimeric,splicing";
$analyze_params[] = "--log ${folder}/analyze_rna_".date('YmdHis',mktime()).".log";
$analyze_params[] = $add_params;

$commands = array("php ".$basedir."Pipelines/analyze_rna.php ".implode(" ", $analyze_params));
$working_directory = realpath($folder."/..");
if ($dryrun) {
	print_r($commands[0]."\n");
}
elseif ($noqueue) {
	chdir($working_directory);
	$parser->execTool("Pipelines/analyze_rna.php", implode(" ", $analyze_params));
}
else {
	$parser->jobsSubmit($commands, $working_directory, get_path('queues_high_mem'), false);
}

?>