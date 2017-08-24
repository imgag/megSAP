<?php

/**
	@page rc_featurecounts
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("rc_featurecounts", "Perform read counting for aligned reads using featureCount contained in the subread package.");
$parser->addInfile("in",  "BAM input file.", false, true);
$parser->addOutfile("out", "Raw count output TSV file.", false);

//optional parameters
$parser->addString("gtf_file", "GTF file containing feature annotations used for read counting.", true, get_path("data_folder")."/dbs/gene_annotations/GRCh37.gtf");
$parser->addString("feature_type", "Feature type used for mapping reads to features.", true, "exon");
$parser->addString("gtf_attribute", "Attribute used as feature ID.", true, "gene_id");

$parser->addEnum("library_type", "Specify the library type, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addFlag("single_end", "Single-end data.");

$parser->addInt("min_mapq", "Minimal mapping quality.", true, 3);

$parser->addFlag("overlap", "Count reads multiple times if they overlap more than one feature.");

$parser->addFlag("keep_summary", "Keep summary file.");

$parser->addInt("threads", "Number of threads used for read counting", true, 4);
extract($parser->parse($argv));

$strandedness = array(
	"unstranded" => 0,
	"reverse" => 2,
	"forward" => 1
);

$tmp_dir = $parser->tempFolder();
$tmp_out = $parser->tempFile();

//build arguments array
$args = array(
	"-a", $gtf_file,
	"-t", $feature_type,
	"-g", $gtf_attribute,
	"-Q", $min_mapq,
	"-T", $threads,
	"-s", $strandedness[$library_type],
	"--tmpDir", $tmp_dir,
	"-o", $tmp_out
);

if(!$single_end) $args[] = "-p -B";
if($overlap) $args[] = "-O";
$args[] = $in;

//execute command
$parser->exec(get_path("feature_counts"), implode(" ", $args), true);

// copy output
$parser->exec("cp", "$tmp_out $out", true);

if ($keep_summary)
{
	$out_summary = dirname($out)."/".basename($out, ".tsv")."_summary.tsv";
	$parser->exec("cp", "{$tmp_out}.summary {$out_summary}", true);
}
