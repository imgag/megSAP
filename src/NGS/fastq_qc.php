<?php

/**
	@page fastq_qc
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("fastq_qc", "Quality control for FASTQ files.");

// supported QC tools
$tools_all = array("readqc","fastqc");

// mandatory arguments
$parser->addString("folder", "Sample folder containing FASTQ files.", false);

// optional arguments
$parser->addString("tools", "Comma-separated list of QC tools to run.", true, implode(",", $tools_all));
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
$parser->addString("out_folder", "Folder where QC results are stored, defaults to sample folder.", true, "default");
$parser->addFlag("import", "Import QC data into NGSD.");
extract($parser->parse($argv));

// check steps
$tools = explode(",", $tools);
$invalid_tools = array_diff($tools, $tools_all);
if (count($invalid_tools) > 0) trigger_error("Unsupported QC tools: ".implode(",",$invalid_tools), E_USER_ERROR);


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
if (in_array("readqc", $tools)) {
	if (count($in_rev) == 0) {
		$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 $in_for_s -out $readqc_out", true);
	} else {
		$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 $in_for_s -in2 $in_rev_s -out $readqc_out", true);
	}

	// if import set and sample id available, import ReadQC results into database
	if ($import && ($name != "")) {
		$parser->execTool("NGS/db_import_qc.php", "-id $name -files $readqc_out -force");
	}
}

// FastQC
if (in_array("fastqc", $tools)) {
	$parser->exec(get_path("fastqc"), "--threads $threads $in_for_s $in_rev_s", true);
}

?>
