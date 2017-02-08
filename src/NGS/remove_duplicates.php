<?php 

/** 
	@page remove_duplicates
*/


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("remove_duplicates", "Removed duplicate reads.");
$parser->addInfile("in",  "Input BAM file.", false);
$parser->addOutfile("out",  "Output BAM file.", false);
extract($parser->parse($argv));

// perform removal of duplicates
$metricsfile = $parser->tempFile("_remove_dup_metrics.txt");
$parser->exec(get_path("picard_duplicates"), "I=$in O=$out M=$metricsfile AS=true", true);

// log metrics file
$file = file($metricsfile);
$parser->log("Metrics file:", $file);

?>
