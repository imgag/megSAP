<?php
/** 
	@page db_background_jobs
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_background_jobs", "Runs jobs in the SGE queue without interferring with the NGSD queuing system.");
$parser->addInfile("commands", "Input text file containing command per line.", false);
$parser->addInt("slots_per_job", "Number of SGE slots to use per command", false);
//optional
$parser->addFloat("max_slots", "Maximum percentage of SGE slots to use.", true, 50.0);
$parser->addInt("sleep_secs", "Number of seconds to sleep between tries to start jobs.", true, 120);
extract($parser->parse($argv));

//init
$user = trim(exec('whoami'));
$queues = get_path("queues_default");

//load commands
$file = file($commands);
$commands = array();
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="")  continue;
	$commands[] = $line;
}

while(count($commands)>0)
{
	print date("Y-m-d h:i:s")."\n";
	
	//determine slots
	list($stdout) = exec2("qstat -u '*' -f | grep default");
	$slots_overall = 0;
	$slots_used = 0;
	foreach($stdout as $line)
	{
		$line = preg_replace('/\s+/', ' ', $line);
		$line = explode(" ", $line);
		list(, $used, $overall) = explode("/", $line[2]);
		$slots_overall += $overall;
		$slots_used += $used;
	}

	$slots_max_used = (int)($max_slots/100.0*$slots_overall);
	$slots_free = $slots_max_used - $slots_used;
	print "  Slots - overall:{$slots_overall} used:{$slots_used} max usable:{$slots_max_used} available:{$slots_free}\n";
	if ($slots_free < $slots_per_job)
	{
		print "  No slots available\n";
	}
	else
	{
		$id = 0;
		while ($slots_free >= $slots_per_job)
		{
			++$id;
			$command = array_shift($commands);
			print "  Starting command: {$command}\n";
			if (contains($command, "megSAP/src/NGS/db_queue_analysis.php")) //queue via NGSD
			{
				exec2($command);
			}
			else //submit to SGE directly
			{
				$sge_folder = get_path("data_folder")."/sge/background_jobs/";
				$base = "{$sge_folder}".date("Ymdhis")."_".str_pad($id, 3, '0', STR_PAD_LEFT)."_{$user}";
				$sge_out = "{$base}.out";
				$sge_err = "{$base}.err";
				$command_sge = "qsub -V -pe smp {$slots_per_job} -b y -wd {$sge_folder} -m n -M ".get_path("queue_email")." -e {$sge_err} -o {$sge_out} -q {$queues} -shell n";
				list($stdout, $stderr) = exec2($command_sge." ".$command);
				
				$sge_id = explode(" ", $stdout[0])[2];
				print "    SGE job id: {$sge_id}\n";
				print "    SGE stdout: {$sge_out}\n";
				print "    SGE stderr: {$sge_err}\n";
			}
			$slots_free -= $slots_per_job;
		}
		print "  Commands remaining: ".count($commands)."\n";
	}
	
	sleep($sleep_secs);
}

?>