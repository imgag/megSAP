<?php

/**
	@page qc_fastq
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("qc_fastq", "Quality control for FASTQ files.");

// mandatory arguments
$parser->addString("folder", "Sample folder containing FASTQ files.", false);

// optional arguments
$parser->addEnum("tool", "QC tool to run.", true, array("readqc","fastqc"), "readqc");
$parser->addInt("threads", "The maximum number of threads to use (for supported tools).", true, 4);
$parser->addString("out_folder", "Folder where QC results are stored, defaults to sample folder.", true, "default");
$parser->addFlag("import", "Import QC data into NGSD.");
extract($parser->parse($argv));

// resolve output folder
if ($out_folder=="default")
{
	$out_folder = $folder;
}

// input FASTQ files
$in_for = glob($out_folder."/*_R1_001.fastq.gz");
$in_rev = glob($out_folder."/*_R2_001.fastq.gz");

if ((count($in_for) == 0) && (count($in_rev) == 0)) {
	trigger_error("No FASTQ files found!", E_USER_ERROR);
}

// output file names
// use sample id, if possible
$p = basename(realpath($out_folder));
$name = "";
if (preg_match('/^Sample_(.+)/', $p, $matches)) {
	$name = $matches[1];
}
$readqc_out  = $out_folder."/".$name."_stats_fastq.qcML";

$in_for_s = implode(" ", $in_for);
$in_rev_s = implode(" ", $in_rev);

// ReadQC
if ($tool == "readqc") {
	if (count($in_rev) == 0) {
		$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 $in_for_s -out $readqc_out", true);
	} else {
		$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 $in_for_s -in2 $in_rev_s -out $readqc_out", true);
	}

	// if import set and sample id available, import ReadQC results into database
	if ($import && ($name != "")) {
		$parser->execTool("NGS/db_import_qc.php", "-id $name -files $readqc_out");
	}
}

// FastQC
if ($tool == "fastqc") {
	$parser->exec(get_path("fastqc"), "--threads $threads $in_for_s $in_rev_s", true);
}

?>