<?php
/** 
	@page export_system 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_system", "Exports a processing system INI file.");
$parser->addString("name", "Processing system short name.", false);
$parser->addOutfile("out", "Output file in INI format.", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//check output file does not exist
if (file_exists($out))
{
	trigger_error("Output file '$out' already exists!", E_USER_ERROR);
}

//export system
$db_conn = DB::getInstance($db);
store_system($db_conn, $name, $out);


?>