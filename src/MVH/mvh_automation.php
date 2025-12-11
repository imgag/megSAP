<?php
/** 
	@page mvh_automation
	
	Queues one export (KDK or GRZ) by adding it to the queue, if the queue has a free slot.
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function running_jobs()
{
	$output = [];
	
	list($stdout) = exec2("qstat -q mvh_srv005", false);
	foreach($stdout as $line)
	{
		$line = trim($line);
		if ($line=="" || starts_with($line, "-----") || starts_with($line, "job-ID")) continue;
		list($id) = explode(" ", $line, 2);
		
		list($stdout2) = exec2("qstat -j {$id} | grep job_args");
		$parts = explode(",", $stdout2[0]);
		foreach($parts as $part)
		{
			if (is_numeric($part)) $output[] = $part;
		}
	}
	
	return $output;	
}

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
}

//parse command line arguments
$parser = new ToolBase("mvh_automation", "Automates MVH upload to GRZ/KDK.");
extract($parser->parse($argv));

//check if queue slots are free
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

//queue
$exclude = running_jobs();
while($available>0)
{
	//check if there is a KDK export  to do (we do KDK exports first, because they are fast)
	$db_mvh = DB::getInstance("MVH");
	$cm_ids = $db_mvh->getValues("SELECT cd.cm_id FROM case_data cd, submission_kdk_se sub WHERE sub.case_id=cd.id AND sub.status='pending'");
	foreach($cm_ids as $cm_id)
	{
		if (in_array($cm_id, $exclude)) continue;
		
		queue_job($queue, $cm_id, "KDK_SE");
		$exclude[] = $cm_id;
		$available -= 1;
		continue(2);
	}
	
	//check if there is a GRZ export to do 
	$cm_ids = $db_mvh->getValues("SELECT cd.cm_id FROM case_data cd, submission_grz sub WHERE sub.case_id=cd.id AND sub.status='pending'");
	foreach($cm_ids as $cm_id)
	{
		if (in_array($cm_id, $exclude)) continue;
		
		queue_job($queue, $cm_id, "GRZ");
		$exclude[] = $cm_id;
		$available -= 1;
		continue(2);
	}
	
	break; //terminate if there are slots available, but no more jobs
}

?>
