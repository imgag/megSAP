<?php
/** 
	@page db_delete_secondary_analysis
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_delete_secondary_analysis", "Deletes a secondary analysis.");
$parser->addString("gsvar", "Substring to match in GSvar file.", true);
$parser->addString("type", "Type to match.", true);
$parser->addFlag("multi", "Allow deletion of several matches.");
$parser->addFlag("debug", "Do not preform deletion, just print commands that would be executed.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

$db = DB::getInstance($db);

//create constrains
$constraints = [];
if ($gsvar!="")
{
	$constraints[] = "gsvar_file LIKE '%{$gsvar}%'";
}
if ($type!="")
{
	$constraints[] = "type='{$type}'";
}
if (count($constraints)==0) trigger_error("Please provide parameters!", E_USER_ERROR);

//get IDs
$ids = $db->getValues("SELECT id FROM secondary_analysis WHERE ".implode(" AND ", $constraints));
$count = count($ids);
if ($count>1 && !$multi) trigger_error("{$count} secondary analysis to delete: ".implode(", ", $ids)." Use '-multi' flag to perform deletion of multiple matches!", E_USER_ERROR);

//delete
foreach($ids as $id)
{
	//folder
	$gsvar = $db->getValue("SELECT gsvar_file FROM secondary_analysis WHERE id={$id}");
	$folder = dirname($gsvar);
	print "Deleting folder $folder\n";
	exec2("rm -rf $folder");
	
	//secondary analysis table entry
	print "Deleting NGSD entry with ID $id\n";
	$db->executeStmt("DELETE FROM secondary_analysis WHERE id={$id}");
}
?>