<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("filter_variants", "Wrapper for VariantFilterAnnotations.");
$parser->addInfile("in", "Input GSvar file.", false);
$parser->addInfile("filters", "Filter definition TXT file.", false);
$parser->addOutfile("out", "Output GSvar file.", false);
$parser->addString("requires_db", "Check if the given DB is enabled before performing the filtering.", true);
extract($parser->parse($argv));

//check for DB
if ($requires_db!="" && !db_is_enabled($requires_db))
{
	trigger_error("Skipped - database '$requires_db' is not enabled!", E_USER_NOTICE);
	exit(0);
}

//exec
exec2(get_path("ngs-bits")."/VariantFilterAnnotations -in {$in} -out {$out} -filters {$filters}");

?>