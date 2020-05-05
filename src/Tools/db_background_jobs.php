<?php
/** 
	@page db_background_jobs
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_background_jobs", "Runs jobs in the SGE queue without interferring with the NGSD queuing system.");
$parser->addInfile("commands", "Input text file containing one 'db_queue_analysis' command per line.", false);
$parser->addInt("slots_per_job", "Number of SGE slots to use per command", false);
//optional
$parser->addFloat("max_slots", "Maximum percentage of SGE slots to use.", true, 80.0);
$parser->addInt("sleep_secs", "Number of seconds to sleep between tries to start jobs.", true, 120);
extract($parser->parse($argv));

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
		while ($slots_free >= $slots_per_job)
		{
			$command = array_shift($commands);
			print "  Starting command: {$command}\n";
			exec2($command);
			$slots_free -= $slots_per_job;
		}
		print "  Commands remaining: ".count($commands)."\n";
	}
	
	sleep($sleep_secs);
}

?>