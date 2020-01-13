<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("compare_variants", "Wrapper for SampleDiff.");
$parser->addInfile("in1", "First GSvar file.", false, false);
$parser->addInfile("in2", "Second GSvar file.", false, false);
$parser->addString("requires_db", "Check if the given DB is enabled before performing the filtering.", true);
extract($parser->parse($argv));

//check for DB
if ($requires_db!="" && !db_is_enabled($requires_db))
{
	trigger_error("Skipped - database '$requires_db' is not enabled!", E_USER_NOTICE);
	exit(0);
}

//exec
list($stdout, $stderr) = exec2(get_path("ngs-bits")."/SampleDiff -window 0 -in1 {$in1} -in2 {$in2} 2>&1");
foreach($stdout as $line)
{
	print $line."\n";
}
foreach($stderr as $line)
{
	print $line."\n";
}

?>