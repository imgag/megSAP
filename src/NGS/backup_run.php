<?php 

/** 
	@page backup_run
*/


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("backup_run", "Creates a backup of a Illumina run folder that was created by RTA 1.9 or higher.");
$parser->addInfile("in",  "Input run folder.", false);
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

//zip archive with full compression
print date("Y-m-d H:i:s")." creating tar.gz file\n";
$tmp_zip = $parser->tempFile(".tar.gz");
$parser->exec("gzip", "-c -9 $tmp_tar > $tmp_zip", true);

//test zip archive integrity
print date("Y-m-d H:i:s")." testing tar.gz file integrity\n";
$parser->exec("gzip", "-t $tmp_zip", true);
$parser->deleteTempFile($tmp_tar);

//calculate MD5 checksum of zip archive
print date("Y-m-d H:i:s")." calulating md5sum checksum of tar.gz file\n";
list($stdout) = $parser->exec("md5sum", "$tmp_zip", true);
list($md5) = explode(" ", $stdout[0]);

//copy zip file to archive folder
print date("Y-m-d H:i:s")." copying tar.gz file into archive\n";
$zipfile = $out_folder."/".basename($in).".tar.gz";
$parser->copyFile($tmp_zip, $zipfile);

//validate remote checksum
print date("Y-m-d H:i:s")." calulating md5sum checksum of copy\n";
list($stdout) = $parser->exec("md5sum", "$zipfile", true);
list($md5_2) = explode(" ", $stdout[0]);
if ($md5 != $md5_2)
{
	trigger_error("MD5 checksums do not match: '$md5' != '$md5_2'!", E_USER_ERROR);
}

//copy log file to remote archive folder
$logfile = $out_folder."/".basename($in).".log";
$parser->copyFile($parser->getLogFile(), $logfile);

//print log file to console
print date("Y-m-d H:i:s")." done\n";
print implode("", file($parser->getLogFile()));

?>
