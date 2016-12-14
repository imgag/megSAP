<?php 
/** 
	@page queue_sample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");


$parser = new ToolBase("queue_sample", "Queues sample for analysis.");
$parser->addString("sample",  "Processed sample identifier.", false);
$parser->addString("steps", "Comma-separated list of processing steps to perform.", true, "ma,vc,an,db,cn");
$parser->addFlag("backup", "Backup files from previous analysis.");
$parser->addFlag("high_priority", "Assign a high priority to job.");
$parser->addString("user", "name of the user who queued the analysis.", true, "");
extract($parser->parse($argv));

//check if the given processed sample is already in queue
function check_ps_in_queue($processed_sample_name)
{
	//build jobID-list for processed sample
	$job_info_file = file(get_path("sample_status_folder")."/jobinfo.txt");
	$ps_name_to_jobIDs=array();
	foreach($job_info_file as $line)//build a row for each entry in jobinfo.txt
	{
		list($jobnumber,,,,$found_proc_sample_name)=explode("\t",$line);
		if ($found_proc_sample_name==$processed_sample_name)
		{
			$ps_name_to_jobIDs[]=$jobnumber;
		}
	}
	
	//check if a job corresponding to the processed sample is in queue
	foreach($ps_name_to_jobIDs as $jobID)
	{
		$qstat_output="";
		$qstat_output=shell_exec("qstat -j ".$jobID." 2> /dev/null"); //don't show error messages
		if ($qstat_output!="")//if the job was found
		{
			trigger_error("Processed Sample ".$processed_sample_name. " is already in queue with jobnumber=".$jobID."!", E_USER_ERROR);
		}
	}
}

//get sample information from NGSD
function get_parameters($processed_sample_name)
{
	//get sample name and process ID
	list($sample_name,$process_id) = explode("_",$processed_sample_name);
	$process_id=ltrim($process_id,'0');
	$db_connect = DB::getInstance('NGSD');
	
	//get tumour status
	$result = $db_connect->executeQuery("SELECT tumor FROM  sample WHERE name='".$sample_name."'");
	$tumor_status=$result[0]['tumor'];	
	
	//get sample id
	$result = $db_connect->executeQuery("SELECT id FROM  sample WHERE name='".$sample_name."'");
	$sample_id=$result[0]["id"];
	
	//get processed sample id and project id
	$result = $db_connect->executeQuery("SELECT id,project_id FROM  processed_sample WHERE sample_id=".$sample_id." AND process_id= ".$process_id);
	if (count($result)==0)
	{
		trigger_error("Could not find processed sample '$processed_sample_name' in NGSD!", E_USER_ERROR);
	}
	$proc_sample_id=$result[0]["id"];
	$projectID=$result[0]["project_id"];
	
	//get project name
	$result = $db_connect->executeQuery("SELECT  name,type FROM project WHERE id=".$projectID);
	$projectname=$result[0]["name"];
	$project_type=$result[0]["type"]; 

	//get run ID
	$result = $db_connect->executeQuery("SELECT sequencing_run_id FROM  processed_sample WHERE sample_id=".$sample_id." AND process_id= ".$process_id);
	$run_id=$result[0]["sequencing_run_id"];

	//get run name
	$result = $db_connect->executeQuery("SELECT name FROM  sequencing_run WHERE id=".$run_id);
	$run_name=$result[0]["name"];
	
	//get WGS
	$result = $db_connect->executeQuery("SELECT sys.type FROM processed_sample ps, processing_system sys WHERE ps.id='{$proc_sample_id}' AND sys.id=ps.processing_system_id");
	$wgs = ($result[0]["type"]=="WGS");
	
	return array($proc_sample_id, $projectID, $process_id, $project_type, $run_name, $run_id, $tumor_status, $projectname, $wgs);
}

//main script starts here!
check_ps_in_queue($sample);

list($proc_sample_id, $projectID, $process_id, $project_type, $run_name, $run_id, $tumor, $project_name, $wgs) = get_parameters($sample);

//determine project folder
$project_folder = get_path("project_folder")."/".$project_type."/".$project_name."/";

//check sample folder is available
$sample_folder = $project_folder."Sample_".$sample."/";
if(!is_dir($sample_folder))
{
	trigger_error("Could not find sample folder '$sample_folder'!", E_USER_ERROR);
}

//check that fastq files are there
$files = glob($sample_folder.$sample."*.fastq.gz");
if(count($files)<2)
{
	trigger_error("Could not find at least two FASTQ files starting with '$sample' in '$sample_folder'!", E_USER_ERROR);
}

//determine command and arguments
if ($tumor)
{
	$outfolder = $project_folder."/Sample_".$sample."/";
	if (!file_exists($outfolder)) mkdir($outfolder);
	
	$command = "php ".repository_basedir()."/src/Pipelines/somatic_dna.php";
	//remove 'cn' step if present
	$steps=str_replace (array(",cn","cn") , "" , $steps);
	$args = "-p_folder {$project_folder} -t_id {$sample} -n_id na -o_folder {$outfolder} -steps {$steps} --log {$outfolder}somatic_dna_".date("Ymdhis").".log";
}
else
{
	$command = "php ".repository_basedir()."/src/Pipelines/analyze.php";
	$args = "-folder {$sample_folder} -name {$sample} -steps {$steps} --log {$sample_folder}analyze_".date("Ymdhis").".log";
	if ($backup)
	{
		$args .= " -backup";
	}
	if ($wgs) //whole genome => use 5 threads (default is 2)
	{
		$args .= " -threads 5";
	}	
}

//submit to queue
$sample_status = get_path("sample_status_folder")."/data/";
$queue = $high_priority ? "-q NGSlong,re_analysis,srv016_long,srv018long" : "-q NGSlong,srv016_long,srv018long";
$slots = $wgs ? "-pe smp 2" : ""; //use two instead of one slot for WGS
$queue_command = "qsub -V $slots -b y -wd $project_folder -m n -M florian.lenz@med.uni-tuebingen.de -e $sample_status $queue -o $sample_status $command $args";
$qsub_return_line = shell_exec($queue_command);

//extract job number
$exploded_qsub_return_line = explode(" ", $qsub_return_line);
$jobnumber = intval($exploded_qsub_return_line[2]);
if($jobnumber>0)
{
	if ($user)
	{
		$user_name=$user;
	}
	else
	{
		$user_name = trim(shell_exec('whoami'));
	}
	$outputline = array($jobnumber,date("d-m-Y_H:i:s"), $run_name, $run_id, $sample, $proc_sample_id, $project_name, $projectID, $user_name);
	$outputline[] = $high_priority ? "high" : "low";
	file_put_contents(get_path("sample_status_folder")."/jobinfo.txt", implode("\t", $outputline)."\n", FILE_APPEND);
}
else
{
	trigger_error("Could not queue job - job number lower than 1 was returned!", E_USER_ERROR);
}