<?php
/** 
	@page mvh_grz_export
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function json_patient()
{
	$output = [];
	return $output;
}

//parse command line arguments
$parser = new ToolBase("mvh_kdk_se_export", "KDK-SE export for Modellvorhaben.");
$parser->addInt("case_id", "'id' in 'data_data' of 'MVH' database.", false);
$parser->addFlag("clear", "Clear export and QC folder before running this script.");
$parser->addFlag("test", "Test mode.");
extract($parser->parse($argv));

//init
$db_mvh = DB::getInstance("MVH");
$db_ngsd = DB::getInstance("NGSD");
$mvh_folder = get_path("mvh_folder");

//check that case ID is valid
$id = $db_mvh->getValue("SELECT id FROM case_data WHERE id='{$case_id}'");
if ($id=="") trigger_error("No case with id '{$case}' in MVH database!", E_USER_ERROR);

//create export folder
$folder = realpath($mvh_folder)."/kdk_se_export/{$case_id}/";
if ($clear) exec2("rm -rf {$folder}");
exec2("mkdir -p {$folder}");
print "export folder: {$folder}\n";

//create JSON
$json = [
	"patient"=>json_patient();
	];

//write JSON
file_put_contents("{$folder}/metadata/metadata.json", json_encode($json, JSON_PRETTY_PRINT));

?>
