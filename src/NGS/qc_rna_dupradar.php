<?php

/**
	@page qc_rna_dupradar
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";


error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("qc_rna_dupradar", "Calculate duplication rate quality control data on RNA-seq sample using dupRadar.");

// mandatory arguments
$parser->addString("folder", "Sample folder containing BAM file.", false);

$parser->addFlag("se", "Flag indicating single-ended sequencing.");
$parser->addInfile("system",  "Processing system INI file (determined from folder name by default).", true);
$parser->addString("out_folder", "Output folder for QC results, defaults to sample-folder/qc/dupradar", true, "");
$parser->addInt("threads", "Number of parallel threads.", true, 4);
extract($parser->parse($argv));

// resolve and create  output folder
$folder = realpath($folder);
if ($out_folder == "") {
	// set default output folder
	$out_folder = $folder."/qc/dupradar/";
}

if (!is_dir($out_folder) && !mkdir($out_folder, 0777, true)) {
	trigger_error("Could not create output folder '".$out_folder."'!", E_USER_ERROR);
}

// extract processed sample id
if (!preg_match('/^Sample_(.+)_(.+)/', basename($folder), $matches)) {
	trigger_error("Could not resolve processed sample from folder name!", E_USER_ERROR);
}
$sample = $matches[1];
$processed_sample = $matches[1]."_".$matches[2];

// extract processing system information
$sys = load_system($system, $processed_sample);
if (!isset($sys['build']) || $sys['build']=="") {
	trigger_error("Build not specified!", E_USER_ERROR);
} else {
	$build = $sys['build'];
}

// dupRadar invocation
$bam = "{$folder}/{$processed_sample}.bam";
$gtf = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

$prog = get_path("dupradar");

// dupRadar parameters
$params = array();
$params[] = $bam;
$params[] = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

// TODO stranded=1 or =2?
if (!isset($sys['stranded']) || (isset($sys['stranded']) && $sys['stranded']>0)) {
	$params[] = "stranded=reverse";
} else {
	$params[] = "stranded=no";
}
if (!$se) {
	$params[] = "paired=yes";
} else {
	$params[] = "paired=no";
}

$params[] = "outdir={$out_folder}";
$params[] = "threads={$threads}";
if (isset($sys['gtfAttribute'])) {
	$params[] = "gtfAttribute=".$sys['gtfAttribute'];
} else {
	$params[] = "gtfAttribute=gene_id";
}
if (isset($sys['featureType'])) {
	$params[] = "gtfFeatureType=".$sys['featureType'];
} else {
	$params[] = "gtfFeatureType=exon";
}

$parser->exec($prog, implode(" ", $params)." > {$out_folder}/dupradar.tsv", true);