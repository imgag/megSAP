<?php 

/** 
	@page backup_run
*/


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("backup_run", "Creates a backup of a run folder.");
$parser->addInfile("in",  "Input run folder.", false);
$parser->addString("when",  "Start time in format '20:15' or 'now'.", true, "now");
$parser->addString("out_folder", "Output folder.", true, "/mnt/storage1/raw_data_archive/runs/");
$parser->addString("sav_folder", "Output folder for SAV data.", true, "/mnt/storage3/raw_data/_sav_archive/");
$parser->addFlag("include_raw_signal", "Backup includes non-basecalled POD5 or FAST5 data.");
$parser->addFlag("test", "Perform a backup in test-mode. No data will be moved on the TSM backup.");
extract($parser->parse($argv));

//check that the correct user is executing the script
$user = exec('whoami');
if ($user!="archive-gs" && !$test) //no check in test-mode
{
	trigger_error("Only user 'archive-gs' can execute this script - use 'sudo -u archive-gs php backup_run.php ...'!", E_USER_ERROR);
}

//do not allow default backup location as output for test-mode
if ($test && ($out_folder == "/mnt/storage1/raw_data_archive/runs/")) trigger_error("ERROR: Do not backup to default output path in test-mode!", E_USER_ERROR);

//strip slashes at the end of folder names
$in = realpath($in);
$in = rtrim($in, "/");
$out_folder = rtrim($out_folder, "/");
$sav_folder = rtrim($sav_folder, "/");

//copy SAV data
$basename = basename($in);
if (file_exists("{$in}/InterOp/"))
{
	$sav_out = "{$sav_folder}/{$basename}/";
	print "Copying SAV data to: {$sav_out}\n";
	$parser->exec("mkdir", "-p {$sav_out}", true);
	$parser->exec("cp", "-R {$in}/InterOp/ {$in}/*.xml {$sav_out}", true);
}

//determine output file names
$archive_basename = $basename;
if (!preg_match("/^([0-9]{2})([0-9]{2})([0-9]{2})_/", $basename, $matches) || !checkdate($matches[2], $matches[3], $matches[1])) 
{
	if (preg_match("/^([0-9]{4})([0-9]{2})([0-9]{2})_/", $basename, $matches) || !checkdate($matches[2], $matches[3], $matches[1]))
	{
		$archive_basename = substr($basename, 2);
		print "Notice: Folder name does start with 4-digit year. Cutting year to 2 digits. Resulting archive name: $archive_basename\n";
	}
	else
	{
		$archive_basename = date("ymd")."_".$basename;
		print "Notice: Folder name does not start with date. Prefixing date: $archive_basename\n";
	}
}

$zipfile = "{$out_folder}/{$archive_basename}.tar.gz";
$logfile = "{$out_folder}/{$archive_basename}.log";

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

// get run id
$split_filename = explode("_", $in);
$run_id = "#".strtr(end($split_filename), array('/' => ''));
if(!preg_match("/^#[0-9]{5}$/", $run_id))
{
	trigger_error("ERROR: Run folder doesn't end with a valid run id!", E_USER_ERROR);
}

//create verified tar archive (uncompressed because otherwise verification is not possible)
print date("Y-m-d H:i:s")." creating tar file\n";
$tmp_tar = $parser->tempFile(".tar");

//decide what to backup depending on run type
$run_parameters_xml = $in."/RunParameters.xml";
if(!file_exists($run_parameters_xml))
{
	trigger_error("ERROR: Required RunParameters.xml file is missing in the run folder!", E_USER_WARNING);
	$is_novaseq_x = false;
}
else
{
	$is_novaseq_x = is_novaseq_x_run($run_parameters_xml);
}
if($is_novaseq_x)
{
	trigger_error("NOTICE: NovaSeq X run detected!", E_USER_NOTICE);
	$analyses = array_filter(glob($in."/Analysis/[0-9]"), 'is_dir');
	if(count($analyses) == 0) trigger_error("ERROR: No analysis found! Please make sure the secondary analysis was performed successfully (at least the demultiplexing).", E_USER_ERROR);
	if(count($analyses) > 1) trigger_error("ERROR: Multiple analysis folders found! Please remove all duplicated secondary analysis, but make sure all necessary fastq files are present in the remaining folder!", E_USER_ERROR);
	//NovaSeq X run => backup fastq/ora files
	$parser->exec("tar", "cfW {$tmp_tar} --exclude='{$basename}/Data' --exclude='{$basename}/Unaligned*' --exclude='{$basename}/Thumbnail_Images' --exclude='{$basename}/Images'"
	." --exclude='*.bam' --exclude='*.bam.bai' --exclude='*.cram' --exclude='*.cram.crai'  --exclude='*.vcf.gz'  --exclude='*.gvcf.gz' -C ".dirname($in)." {$basename}", true);
}
else
{
	trigger_error("NOTICE: Normal run detected!", E_USER_NOTICE);
	if ($include_raw_signal)
	{
		$exclude_raw = "";
	}
	else
	{
		$exclude_raw = "--exclude='*.fast5' --exclude='*.pod5'";
	}
	//normal run => backup BCL files
	$parser->exec("tar", "cfW {$tmp_tar} {$exclude_raw} --exclude='{$basename}/Unaligned*' --exclude='{$basename}/Thumbnail_Images' --exclude='{$basename}/Images' -C ".dirname($in)." {$basename}", true);
}


//zip archive with low compression (way faster than full compression and the contents are already compressed in most cases)
print date("Y-m-d H:i:s")." creating tar.gz file\n";
$parser->exec("pigz", "-p 8 -c -1 $tmp_tar > $zipfile", true);

//log file size
print date("Y-m-d H:i:s")." get compressed file size\n";
$parser->exec("du", "-h $zipfile", true);

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

// update NGSD run information
$db_name = ($test)?"NGSD_TEST":"NGSD";
if (db_is_enabled($db_name))
{
	// get run entry from NGSD
	$db = DB::getInstance($db_name);
	$query = "SELECT id FROM sequencing_run WHERE name = '{$run_id}'";
	$res = $db->executeQuery($query);

	if (count($res) == 0)
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
	trigger_error("No Connection to the NGSD available! Could not update run information.", E_USER_WARNING);
}

//print log file to console
print date("Y-m-d H:i:s")." done\n";
print implode("", file($parser->getLogFile()));

?>
