<?php 

/** 
	@page backup_project
*/


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("backup_project", "Creates a backup of a project folder.");
$parser->addInfile("in",  "Input project folder.", false);
//optional
$parser->addString("when",  "Start time in format '20:15' or 'now'.", true, "now");
$parser->addString("out_folder", "Output folder path.", true, "/mnt/storage1/raw_data_archive/projects/");
$parser->addString("tmp_folder", "Temporary folder. If unset a user-specific folder in ".sys_get_temp_dir()." is used.", true, "");
extract($parser->parse($argv));

//check that the correct user is executing the script
$user = exec('whoami');
if ($user!="archive-gs")
{
	trigger_error("Only user 'archive-gs' can execute this script - use 'sudo -u archive-gs php backup_project.php ...'!", E_USER_ERROR);
}

//check that tmp folder is writable for archive-gs
if ($tmp_folder!="" && !is_writable2($tmp_folder."/test"))
{
	trigger_error("Temporary folder '{$tmp_folder}' not writable for 'archive-gs'!", E_USER_ERROR);
}

//strip slashes at the end of folder names
$in = realpath($in);
$in = rtrim($in, "/");
$out_folder = rtrim($out_folder, "/");

//determine output file names
$basename = basename($in);
if (!preg_match("/^([0-9]{2})([0-9]{2})([0-9]{2})_/", $basename, $matches) || !checkdate($matches[2], $matches[3], $matches[1]))
{
	$basename = date("ymd")."_".$basename;
	print "Notice: Folder name does not start with date. Prefixing date: $basename\n";
}
$zipfile = "{$out_folder}/{$basename}.tar.gz";
$logfile = "{$out_folder}/{$basename}.log";

//set temporary logfile if unset
if ($parser->getLogFile()==null)
{
	$log_file = $parser->tempFile(".log", null, $tmp_folder);
	$parser->setLogFile($log_file);
	print "Log file: $log_file\n";
}

//wait until start time
if ($when!="now")
{
	$ts = strtotime($when);
	if ($ts < time()) $ts += 24*60*60; //plus one day
	$sleep = $ts - time();
	print "Waiting until $when (".time_readable($sleep).")\n";
	sleep($sleep);
}

//create verified tar archive (uncompressed because otherwise verification is not possible)
$tmp_tar = $parser->tempFile(".tar", null, $tmp_folder);
print "Temporary tar file: $tmp_tar\n";
print date("Y-m-d H:i:s")." creating tar file\n";
$parser->exec("tar", "cfW $tmp_tar -C ".dirname($in)." ".basename($in), true);

//zip archive with low compression (way faster than full compression and the contents are already compressed in most cases)
print date("Y-m-d H:i:s")." creating tar.gz file\n";
$parser->exec("pigz", "-p 8 -c -1 $tmp_tar > $zipfile", true);

//test zip archive integrity
print date("Y-m-d H:i:s")." testing tar.gz file integrity\n";
$parser->exec("pigz", "-p 8 -t $zipfile", true);
$parser->deleteTempFile($tmp_tar);

//calculate MD5 checksum of zip archive
print date("Y-m-d H:i:s")." calulating md5sum checksum of tar.gz file\n";
list($stdout) = $parser->exec("md5sum", "$zipfile", true);
list($md5) = explode(" ", $stdout[0]);

//copy log file to remote archive folder
$parser->copyFile($parser->getLogFile(), $logfile);

//print log file to console
print date("Y-m-d H:i:s")." done\n";
print implode("", file($parser->getLogFile()));

//flag project as archived if in NGSD
$project_name = basename($in);
$db = DB::getInstance("NGSD");
$id = $db->getValue("SELECT id FROM project WHERE name='{$project_name}'", -1);
if ($id!=-1)
{
	$db->executeStmt("UPDATE project SET archived='1' WHERE name='{$project_name}'");
	print date("Y-m-d H:i:s")." flagged project as archived in NGSD!\n";
}

?>
