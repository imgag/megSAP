<?php

/** 
	@page db_init
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_init", "Initializes the NGSD.");
$parser->addFlag("test", "Work on test database instead of production database");
//optional
$parser->addInfile("add", "Text file with SQL queries to import at the end of the initialization.", true);
extract($parser->parse($argv));


$db = DB::getInstance($test ? "NGSD_TEST" : "NGSD");

print "Initializing NGSD...\n";

$tables = $db->getValues("SHOW TABLES");
if (count($tables)>0)
{
	print "  skipped because NGSD is already initialized: ".count($tables)." tables found!\n";
}
else
{
	$args = [];
	$in_files = [];
	if ($test===true)
	{
		$args[] = "-test";
	}
	if ($add!="")
	{
		$args[] = "-add {$add}";
		$in_files[] = $add;
	}
	$parser->execApptainer("ngs-bits", "NGSDInit", implode(" ", $args), $in_files);
	print "  done\n";
}

?>
