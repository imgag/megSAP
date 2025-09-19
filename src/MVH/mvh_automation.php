<?php
/** 
	@page mvh_automation
	
	Queues one export (KDK or GRZ) by adding it to the queue, if the queue has a free slot.
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);


function queue_job($queue, $cm_id, $type)
{
	print "queuing export of '$cm_id' to '$type'\n";
	
	//get email
	$db = DB::getInstance("NGSD");
	$email = $db->getValue("SELECT email FROM user WHERE user_id='ahsturm1'");
	
	//create command
	$sge_folder = get_path("data_folder")."/sge/mvh/";
	$base = "{$sge_folder}".date("Ymdhis")."_{$cm_id}_{$type}";
	$sge_out = "{$base}.out";
	$sge_err = "{$base}.err";
	$command_sge = "qsub -V -pe smp 1 -b y -wd {$sge_folder} -m a -M {$email} -e {$sge_err} -o {$sge_out} -q {$queue}";
	$command = "{$command_sge} php ".repository_basedir()."/src/MVH/mvh_wrapper.php -cm_id {$cm_id} -type {$type}";
	print "  command: {$command}\n";
	print "  SGE STDERR file: {$sge_err}\n";
	print "  SGE STDOUT file: {$sge_out}\n";
	
	//queue command
	list($stdout, $stderr) = exec2($command);
	
	//print JOB infos
	foreach($stdout as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		print "  STDOUT: {$line}\n";
	}
	foreach($stderr as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		print "  STDERR: {$line}\n";
	}
	
	//terminate execution (we only submit one job)
	exit(0);
}

//parse command line arguments
$parser = new ToolBase("mvh_automation", "Automates MVH upload to GRZ/KDK.");
extract($parser->parse($argv));


//check if queue is empty, abort otherwise
list($stdout) = exec2("qstat -g c | grep mvh_");
if (count($stdout)!=1) trigger_error("Could not determine MVH queue status!", E_USER_ERROR);
list($queue, $load, $used, $res, $available, $total) = explode(" ", trim(preg_replace('/\s+/', ' ', $stdout[0])));
print "queue: {$queue}\n"; 
print "queue slots available: {$available}\n"; 
if ($available==0)
{
	print "Queue '$queue' is busy or deactivated > nothing to do\n";
	exit(0);
}

//check if there is a KDK export  to do (we do KDK exports first, because they are fast)
$db_mvh = DB::getInstance("MVH");
$cm_ids = $db_mvh->getValues("SELECT cd.cm_id FROM case_data cd, submission_kdk_se sub WHERE sub.case_id=cd.id AND sub.status='pending'");
foreach($cm_ids as $cm_id)
{
	queue_job($queue, $cm_id, "KDK_SE");
}
//check if there is a GRZ export to do 
$cm_ids = $db_mvh->getValues("SELECT cd.cm_id FROM case_data cd, submission_grz sub WHERE sub.case_id=cd.id AND sub.status='pending'");
foreach($cm_ids as $cm_id)
{
	queue_job($queue, $cm_id, "GRZ");
}

?>
