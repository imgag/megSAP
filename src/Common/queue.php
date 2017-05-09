<?php 

require_once("genomics.php");

function serious_error($error_file)
{
	$error_file_lines=file($error_file);
	foreach($error_file_lines as $line)
	{
		if ((trim($line)!="No entry for terminal type \"unknown\";")&&(trim($line)!="using dumb terminal settings.")) return true;
	}
	return false;
}
	
function jobStatus($jobID)
{
	$status = null;
	
	//(1) check qstat for running jobs
	if(is_null($status))
	{
		exec("qstat -u \"*\" 2>&1",$qstat_output);	//get status of all running jobs of all users from sungrid engine
//		array_shift($qstat_output);	//remove header
//		array_shift($qstat_output);	//remove header
		$job_to_qstate=array();	//build map of jobnumber->status
		for($i=2;$i<count($qstat_output);++$i)
		{
			$qstat_line = $qstat_output[$i];
			$exploded_qstat_line=array_filter(explode(" ",$qstat_line));	//filter since there may be multiple whitespaces between values resulting in multiple empty  fields
			$exploded_qstat_line=array_values($exploded_qstat_line);	//recalculate indexes after filtering
			list($jobnumber_key,,,,$state)=$exploded_qstat_line;
			$job_to_qstate[$jobnumber_key]=$state;
		}

		if (key_exists($jobID, $job_to_qstate))	//looks like job is currently in queue, get status
		{
			$qstate = $job_to_qstate[$jobID];
			switch ($qstate)
			{
				case 'Eqw':
					exec("qstat -j $jobID 2>&1",$qstat_output);
					$status = "queue error: ".implode("\n",$qstat_output);	//error while trying to put analysis in queue
					break;
				case 'qw':
					$status = "queued";	//analysis in queue, but not started
					break;
				case 'r':	
					$status = "running";	//analysis running
					break;
				case 't':
					$status = "transferring";	//analysis was started but temporary suspended due to high server load
					break;
				case 'T':
					$status = "suspended";	//analysis was started but temporary suspended due to high server load (threshold)
					break;
				case 's':
					$status = "suspended";	//analysis was started but temporary suspended due to high server load (qmod)
					break;
				case 'S':
					$status = "suspended";	//analysis was started but temporary suspended due to high server load (queue)
					break;
				default:
					trigger_error("Unknown queue status (".$qstate.").",E_USER_ERROR);
					break;
			}
		}
	}
	
	//(2) check qacct for finished jobs
	if(is_null($status))
	{
		$count = 0;
		$sleep = 0;
		$max_count = 6;
		$min_sleep = 5;
		$exit_status = null;
		while (is_null($exit_status) && $count < $max_count)	//repeat several times since it may take some time until qacct will answer correctly (flush_time in conf master is currently set to 15 seconds)
		{
			if($count>0)
			{
				$sleep += $min_sleep;
				trigger_error("Job-ID: $jobID. Waiting for $sleep seconds...",E_USER_NOTICE);
			}
			sleep($sleep);
						
			$extra = "";
			$acc = get_path("queue_accounting_file");
			if(is_file($acc))	$extra = "-f $acc";
			exec("qacct $extra -j $jobID 2>&1", $qacct_output);
			foreach($qacct_output as $qacct_line)
			{
				if($qacct_line[0]=="=")	continue;	//skip header line

				$exploded_qacct_line=array_filter(explode(" ",$qacct_line), function($item) {return $item!="";});	//filter since there may be multiple whitespaces between values resulting in multiple empty  fields
				$exploded_qacct_line=array_values($exploded_qacct_line);	//recalculate indexes after filtering
				list($key,$value)=$exploded_qacct_line;
				if($key=="jobnumber" && $value!=$jobID)	trigger_error("Found wrong job ID: looked for '$jobID' found '$value'.",E_USER_ERROR);
				if($key=="exit_status")	$exit_status = $value;
			}
			if(!is_null($exit_status))
			{
				if($exit_status==0)	$status = "finished";
				else if($exit_status>0)	$status = "job error";
			}

			++$count;
		}
	}

	//(3) error if unidentified status
	if(is_null($status))
	{
		$qstat = "";
		foreach($job_to_qstate as $j => $q)
		{
			$qstat .= "$j => $q, ";
		}
		trigger_error("Could not determine status for job-ID $jobID. Qstat returned '".(count($qstat_output)==0?"null":implode("\n",$qstat_output))."'; qacct returned '".(count($qacct_output)==0?"null":implode("\n",$qacct_output))."'.",E_USER_ERROR);
	}
	
	return $status;
}

function jobDelete($jobnumber)
{
	return(exec('qdel '.$jobnumber));
}

function jobsWait($jobnumbers)
{
	sleep(5);	//wait a moment because it might take a little time for job to enter queue
	
	$qstats = array();	//collect qstat messages
	
	foreach($jobnumbers as $jobnumber)
	{
		$qstat = jobStatus($jobnumber);
		while(($qstat=="queued")||($qstat=="suspended")||($qstat=="running") || ($qstat=="transferring"))	//sleep while job waiting or running
		{
			sleep(20);
			$qstat = jobStatus($jobnumber);
		}
		$qstats[$jobnumber] = $qstat;		
	}
	
	return $qstats;
}

function jobSubmit($command, $working_directory, $queue, $status_folder, $wait=false)
{
	//build qsub command
	$exec_line="qsub -V -b y -wd $working_directory -m n -M ".get_path("queue_email")." -q $queue -e $status_folder -o $status_folder $command";
	//submit to queue
	$qsub_return_line = system($exec_line);

	//extract job number
	$exploded_qsub_return_line=explode(" ",$qsub_return_line);
	if(!is_numeric($exploded_qsub_return_line[2]))	trigger_error("Could not extract job ID. Return message: ".$qsub_return_line,E_USER_ERROR);
	$jobnumber=intval($exploded_qsub_return_line[2]);

	if(is_numeric($jobnumber))
	{
		$user_name = exec('whoami');
		$command_short = strlen($command)<60 ? $command : substr($command, 0, 60)."...";
		$outputline = array($jobnumber,date("d-m-Y_H:i:s"), $command_short, $working_directory, $user_name."\n");
		file_put_contents(get_path("sample_status_folder")."/jobinfo.txt", implode("\t", $outputline), FILE_APPEND);

		if ($wait)
		{
			jobsWait(array($jobnumber));
		}
	}
	
	return $jobnumber;
}

?>