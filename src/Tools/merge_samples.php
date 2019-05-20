<?php
/** 
	@page merge_samples 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("merge_samples", "Merges a processed sample into another processed sample.");
$parser->addString("ps", "Processed sample to merged into the second sample.", false);
$parser->addString("into", "Processed sample into which the first sample is merged.", false);
extract($parser->parse($argv));

//get processed sample infos
$db_conn = DB::getInstance("NGSD");
$info1 = get_processed_sample_info($ps);
$info2 = get_processed_sample_info($into);

//check input folders
$folder1 = $info1['ps_folder'];
if (!is_dir($folder1))
{
	trigger_error("Sample folder '$folder1' does not exist!", E_USER_ERROR);
}
$folder2 = $info2['ps_folder'];
if (!is_dir($folder2))
{
	trigger_error("Sample folder '$folder2' does not exist!", E_USER_ERROR);
}

//make sure the "+merged_samples" folder exists
$backup_folder = dirname($folder1)."/+merged_samples";
if (!file_exists($backup_folder))
{
	if (!mkdir($backup_folder))
	{
		trigger_error("Could not create backup folder '$backup_folder'!", E_USER_ERROR);
	}
}

//move FASTQ files
exec2("mv $folder1/*.fastq.gz $folder2/");

//move sample folder to backup
exec2("mv $folder1 $backup_folder/");

//NGSD: remove variants/qc for 'ps'
$ps_id = $info1['ps_id'];
$db->executeStmt("DELETE FROM detected_variant WHERE processed_sample_id='$ps_id'");
$db->executeStmt("DELETE FROM processed_sample_qc WHERE processed_sample_id='$ps_id'");

//NGSD: mark samples as merged
$db->executeStmt("INSERT INTO `merged_processed_samples` (`processed_sample_id`, `merged_into`) VALUES ({$ps_id},".$info2['ps_id'].")");

?>