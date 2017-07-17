<?php
/** 
	@page merge_samples 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("merge_samples", "Merges a processed sample into another processed sample.");
$parser->addString("s1", "Sample folder that is merged into the second sample.", false);
$parser->addString("s2", "Sample folder into which the first sample is merged.", false);
extract($parser->parse($argv));

//check input folders
if (!is_dir($s1))
{
	trigger_error("Sample folder '$s1' does not exist!", E_USER_ERROR);
}
if (!is_dir($s2))
{
	trigger_error("Sample folder '$s2' does not exist!", E_USER_ERROR);
}

//make sure the "+merged_samples" folder exists
$backup_folder = dirname($s1)."/+merged_samples";
if (!file_exists($backup_folder))
{
	if (!mkdir($backup_folder))
	{
		trigger_error("Could not create backup folder '$backup_folder'!", E_USER_ERROR);
	}
}

//move FASTQ files
exec2("mv $s1/*.fastq.gz $s2/");

//move sample folder to backup
exec2("mv $s1 $backup_folder/");

//remove detected variants for s1
$name = strtr($s1, array("Sample_"=>""));
$ps_id = get_processed_sample_id($name, false);
if ($ps_id==-1)
{
	print "Skipping NGSD cleanup because sample '$name' was not found!\n";
}
else
{
	$db = DB::getInstance("NGSD");
	$db->executeStmt("DELETE FROM detected_variant WHERE processed_sample_id='$ps_id'");
	$db->executeStmt("DELETE FROM processed_sample_qc WHERE processed_sample_id='$ps_id'");
}


?>