<?php 

/** 
	@page backup_user
*/


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("backup_user", "Creates a backup of a user folder.");
$parser->addInfile("in",  "Input user folder.", false);
$parser->addString("out_folder", "Output folder path.", true, "/mnt/storage1/raw_data_archive/users/");
extract($parser->parse($argv));

//check that the correct user is executing the script
$user = exec('whoami');
if ($user!="archive-gs")
{
	trigger_error("Only user 'archive-gs' can execute this script - use 'sudo -u archive-gs php backup_user.php ...'!", E_USER_ERROR);
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
	$log_file = $parser->tempFile(".log");
	$parser->setLogFile($log_file);
	print "Log file: $log_file\n";
}

//create verified tar archive (uncompressed because otherwise verification is not possible)
print date("Y-m-d H:i:s")." creating tar file\n";
$tmp_tar = $parser->tempFile(".tar");
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

?>
