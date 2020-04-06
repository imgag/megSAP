<?php 

/** 
	@page backup_run
*/


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("backup_run", "Creates a backup of a (sequencing run) folder.");
$parser->addInfile("in",  "Input folder.", false);
$parser->addString("when",  "Start time in format '20:15' or 'now'.", false);
$parser->addString("out_folder", "Output folder path.", false);
extract($parser->parse($argv));

//check that the correct user is executing the script
$user = exec('whoami');
if ($user!="archive-gs")
{
	trigger_error("Only user 'archive-gs' can execute this script!", E_USER_ERROR);
}

//strip slashes at the end of folder names
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
	$log_file = $parser->tempFile(".log");
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
print date("Y-m-d H:i:s")." creating tar file\n";
$tmp_tar = $parser->tempFile(".tar");
$parser->exec("tar", "cfW $tmp_tar --exclude '$in/Data/Intensities/L00?/C*' --exclude '$in/Unaligned*' --exclude '$in/Thumbnail_Images' --exclude '$in/Images' $in/", true);

//zip archive with medium compression (way faster than full compression and the contents are already compressed in most cases)
print date("Y-m-d H:i:s")." creating tar.gz file\n";
$parser->exec("gzip", "-c -5 $tmp_tar > $zipfile", true);

//test zip archive integrity
print date("Y-m-d H:i:s")." testing tar.gz file integrity\n";
$parser->exec("gzip", "-t $zipfile", true);
$parser->deleteTempFile($tmp_tar);

//calculate MD5 checksum of zip archive
print date("Y-m-d H:i:s")." calulating md5sum checksum of tar.gz file\n";
list($stdout) = $parser->exec("md5sum", "$zipfile", true);
list($md5) = explode(" ", $stdout[0]);

//copy log file to remote archive folder
$parser->copyFile($parser->getLogFile(), $logfile);


// update NGSD run information
if (db_is_enabled("NGSD"))
{
	// get run id
	$split_filename = explode("_", $in);
	$run_id = "#".strtr(end($split_filename), array('/' => ''));

	// get run entry from NGSD
	$db = DB::getInstance("NGSD");
	$query = "SELECT id FROM sequencing_run WHERE name = '{$run_id}'";
	$res = $db->executeQuery($query);

	if (sizeof($res) == 0)
	{
		trigger_error("No run with name \"$run_id\" found in NGSD", E_USER_WARNING);
	}
	else
	{
		$query = "UPDATE sequencing_run SET backup_done = 1 WHERE id = '".$res[0]["id"]."'";
		$db->executeStmt($query);
	}
}
else
{
	trigger_error("No Connection to the NGSD available! Couldnt update run information.", E_USER_WARNING);
}

//print log file to console
print date("Y-m-d H:i:s")." done\n";
print implode("", file($parser->getLogFile()));

?>
