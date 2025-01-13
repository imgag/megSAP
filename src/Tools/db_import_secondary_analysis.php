<?php 
/** 
	@page db_import_secondary_analysis
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_secondary_analysis", "Imports secondary analysis into NGSD.");
$parser->addInfile("gsvar",  "Output GSvar file.", false);
$parser->addString("type", "Secondary analysis type.", false);
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//database connection
$db = DB::getInstance($db);

//check type
$valid_types = $db->getEnum("secondary_analysis", "type");
if (!in_array($type, $valid_types))
{
	trigger_error("Invalid analysis type '{$type}'! Valid types are: '".implode("','", $valid_types)."'!", E_USER_ERROR);
}

//normalize and check GSvar
$gsvar = realpath($gsvar);
if (!file_exists($gsvar))
{
	trigger_error("GSvar file '{$type}' does not exist!", E_USER_ERROR);
}

//check analysis does not exist
$id = $db->getValue("SELECT id FROM secondary_analysis WHERE gsvar_file='{$gsvar}'", -1);
if ($id!=-1)
{
	trigger_error("Secondary analysis already exists with id '$id' > Skipping import!", E_USER_NOTICE);
	return;
}

//import
$db->executeStmt("INSERT INTO `secondary_analysis`(`type`, `gsvar_file`) VALUES ('{$type}','{$gsvar}')");

?>
