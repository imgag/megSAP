<?php

/** 
	@page filter_kraken2
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("filter_kraken2", "Wrapper for kraken2.");
$parser->addInfileArray("in", "Input filename(s).", false, true);

//optional
$parser->addOutfile("unclassified_out", "Print unclassified sequences to filename", true);
$parser->addOutfile("classfied_out", "Print classified sequences to filename", true);
$parser->addOutfile("output", "Print output to filename (default: stdout)", true);
$parser->addOutfile("report", "Print a report with aggregate counts/clade to file", true);

$parser->addString("db", "Path to Kraken2 database to use", true, get_path("data_folder")."/dbs/kraken2_filter_hb");
$parser->addInt("threads", "Number of threads to use", true, 1);
$parser->addInt("min_bq", "Minimum base quality used in classification", true, 0);
$parser->addInt("min_hg", "Minimum number of hit groups needed to make a call", true, 2);
$parser->addFloat("confidence", "Confidence score threshold; must be in [0, 1]", true, 0.0);
$parser->addFlag("quick", "Quick operation (use first hit or hits)");
$parser->addFlag("memory_mapping", "Avoids loading database into RAM");
$parser->addFlag("paired", "The filenames provided have paired-end reads");
$parser->addFlag("use_names", "Print scientific names instead of just taxids");
$parser->addFlag("gzip_compressed", "Input files are compressed with gzip");
$parser->addFlag("bzip2_compressed", "Input files are compressed with bzip2");
extract($parser->parse($argv));

//create args array for kraken2 container
$args = [];
$args[] = "--db {$db}";
$args[] = "--threads {$threads}";
$args[] = "--minimum-base-quality {$min_bq}";
$args[] = "--minimum-hit-groups {$min_hg}";
$args[] = "--confidence {$confidence}";

if($quick) $args[] = "--quick";
if($memory_mapping) $args[] = "--memory-mapping";
if($paired) $args[] = "--paired";
if($use_names) $args[] = "--use-names";
if($gzip_compressed) $args[] = "--gzip-compressed";
if($bzip2_compressed) $args[] = "--bzip2-compressed";

$out_files = array();
if ($unclassified_out!="") 
{
	$args[] = "--unclassified-out {$unclassified_out}";
	$out_files[] = $unclassified_out;
}
if ($classfied_out!="") 
{
	$args[] = "--classified-out {$classfied_out}";
	$out_files[] = $classfied_out;
}
if ($output!="")
{
	$args[] = "--output {$output}";
	if ($output!="-") $out_files[] = $output;
}
if ($report!="") 
{
	$args[] = "--report {$report}";
	$out_files[] = $report;
}

// infiles always as the last argument
$args[] = implode(" ", $in);

$in_files = [$db];
$in_files = array_merge($in_files, $in);

//exec kraken2 container
$parser->execApptainer("kraken2", "kraken2", implode(" ", $args), $in_files, $out_files);

?>