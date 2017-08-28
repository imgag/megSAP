<?php

/**
	@page qc_rna_qorts
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";


error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("qc_rna_qorts", "Calculate extensive quality control data on RNA-seq sample using QoRTs.");

// mandatory arguments
$parser->addString("folder", "Sample folder containing BAM file.", false);

$parser->addFlag("se", "Flag indicating single-ended sequencing.");
$parser->addInfile("system",  "Processing system INI file (determined from folder name by default).", true);
$parser->addString("out_folder", "Output folder for QC results, defaults to sample-folder/qc/qorts", true, "");
$parser->addFlag("sort_by_name", "Sort input BAM file by name.");
extract($parser->parse($argv));

// resolve and create  output folder
$folder = realpath($folder);
if ($out_folder == "") {
	// set default output folder
	$out_folder = $folder."/qc/qorts/";
}
else if (substr($out_folder,-1) != "/") {
	// ensure trailing slash (needed by qorts)
	$out_folder = $out_folder."/";
}
if (!is_dir($out_folder) && !mkdir($out_folder, 0777, true)) {
	trigger_error("Could not create output folder!", E_USER_ERROR);
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

// extract read count and raw/original read lengths from NGSD
$db = DB::getInstance("NGSD");
$res = $db->executeQuery("SELECT id, name, name_external FROM sample WHERE name='{$sample}'");
if (count($res) != 1) {
	trigger_error("Sample not found in database or more than one hit!", E_USER_ERROR);
}
list($ngsd_sid, $ngsd_name, $ngsd_name_e) = array_values($res[0]);

$process_id = intval($matches[2]);
$res = $db->executeQuery("SELECT ps.id FROM processed_sample ps WHERE sample_id='$ngsd_sid' AND process_id={$process_id}");

$ngsd_psid = $res[0]['id'];

// number of reads
$res = $db->executeQuery("SELECT qc.value FROM processed_sample_qc qc WHERE processed_sample_id='$ngsd_psid' AND qc_terms_id=16");
if (count($res) != 1) {
	trigger_error("Numer of reads not found in NGSD, QC results will be incomplete.", E_USER_WARNING);
	$read_count = -1;
} else {
	$read_count = $res[0]['value'];
	if (!$se) {
		$read_count = $read_count/2;
	}
}
// raw read length
$res = $db->executeQuery("SELECT qc.value FROM processed_sample_qc qc WHERE processed_sample_id='$ngsd_psid' AND qc_terms_id=28");
if (count($res) != 1) {
	trigger_error("Read length not found in NGSD, QC results will be incomplete.", E_USER_WARNING);
	$read_length = -1;
} else {
	$read_length = $res[0]['value'];
}

// QoRTs invocation
$bam = "{$folder}/{$processed_sample}.bam";
$gtf = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";
$genome = get_path("data_folder")."/genomes/{$build}.fa";

$prog = get_path("qorts");

// QoRTs parameters
$params = array();
$params[] = "QC";

$params[] = "--title {$processed_sample} --outfilePrefix {$processed_sample}_";

// TODO stranded=1 or =2?
if (!isset($sys['stranded']) || (isset($sys['stranded']) && $sys['stranded']>0)) {
	$params[] = "--stranded";
}
if ($se) {
	$params[] = "--singleEnded";
}
if ($read_count > 0) {
	$params[] = "--seqReadCt {$read_count}";
}
if ($read_length > 0) {
	$params[] = "--maxReadLength {$read_length}";
}

if ($sort_by_name) {
	$bam_namesorted = $parser->tempFile(".bam");
	$sort_tmp = $parser->tempFile();
	$parser->exec(get_path("samtools"), "sort -T {$sort_tmp} -n -o {$bam_namesorted} -@ 2 -m 4G {$bam}", true);
	$bam = $bam_namesorted;
}

// $params[] = "--genomeFA {$genome}";

$params[] = "--addFunctions writeGeneBody,makeJunctionBed,calcDetailedGeneCounts";
$params[] = "--skipFunctions NVC,writeClippedNVC";
$params[] = "--generatePlots";

// following parameters must be the last three
$params[] = $bam;
$params[] = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";
$params[] = $out_folder;

$parser->exec($prog, implode(" ", $params), true);
