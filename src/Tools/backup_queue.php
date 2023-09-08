<?php 

/** 
	@page backup_queue
*/


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("backup_queue", "Queues the backup of a folder in SGE.");
$parser->addInfile("in",  "Absolute path to input folder.", false);
$parser->addEnum("mode",  "Mode.", false, array("run", "project", "user"));
$parser->addString("email", "Email used for notification when SGE job has finished (NGSD login also works).", false);
$parser->addFlag("include_raw_signal", "Backup includes non-basecalled POD5 or FAST5 data.");
extract($parser->parse($argv));

//check that the correct user is executing the script
$user = exec('whoami');
if ($user!="archive-gs")
{
	trigger_error("Only user 'archive-gs' can execute this script - use 'sudo -u archive-gs php backup_queue.php ...'!", E_USER_ERROR);
}

//check input folder exists
$in = realpath($in);
if (!file_exists($in))
{
	trigger_error("Input folder '$in' does not exist!", E_USER_ERROR);
}

//check backup script exists
$backup_script = repository_basedir()."/src/Tools/backup_{$mode}.php";
if (!file_exists($backup_script))
{
	trigger_error("Backup script '$backup_script' does not exist!", E_USER_ERROR);
}

//get email from NGSD if necessary
if (!contains($email, "@"))
{
	$db = DB::getInstance("NGSD");
	$email = $db->getValue("SELECT email FROM user WHERE user_id='{$email}'");
}

//create command
$sge_folder = get_path("data_folder")."/sge/archive/";
$base = "{$sge_folder}".date("Ymdhis")."_".basename($in)."_{$user}";
$sge_out = "{$base}.out";
$sge_err = "{$base}.err";
$command_sge = "qsub -V -pe smp 1 -b y -wd {$sge_folder} -m ea -M {$email} -e {$sge_err} -o {$sge_out} -q archive_srv018";
$command = "{$command_sge} php {$backup_script} -in {$in}";
if ($include_raw_signal)
{
	$command = "{$command} -include_raw_signal";
}

print "    SGE command: {$command}\n";

//exucute command
list($stdout, $stderr) = exec2($command);
$sge_id = explode(" ", $stdout[0])[2];
print "    SGE job id: {$sge_id}\n";
print "    SGE stdout: {$sge_out}\n";
print "    SGE stderr: {$sge_err}\n";

?>
