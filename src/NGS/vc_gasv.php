<?php 
/** 
	@page gasv
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("gasv", "Run GASV. Only WGS.");
$parser->addInfile("bam_file",  "Bam file.", false);
$parser->addString("out", "Output folder.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh37");
extract($parser->parse($argv));

chdir($out);	//move to out_folder
$parser->setLogFile(realpath($parser->getLogFile()));				//redirect log file
//TODO fix version numbers for tool execution
//prepare
$bam2gasv = basename($bam_file);
$parser->exec(get_path("GASV_bam2gasv"), "$bam_file -OUTPUT_PREFIX $bam2gasv",true);

//RUN
$bam2gasv .= ".gasv.in";
$parser->exec(get_path("GASV_gasv"), "--batch $bam2gasv",true);