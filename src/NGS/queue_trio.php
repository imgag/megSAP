<?php 
/** 
	@page queue_trio
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("queue_trio", "Queues a trio for analysis.");
$parser->addString("f",  "Processed sample identifier of father", false);
$parser->addString("m",  "Processed sample identifier of mother", false);
$parser->addString("c",  "Processed sample identifier of child", false);
$parser->addFlag("high_priority", "Assign a high priority to job.");
extract($parser->parse($argv));

//check if a trio analysis for that family is already in queue
function check_trio_in_queue($f,$m,$c)
{
	
	//build list of jobs containing family
	$trio=$f."|".$m."|".$c;
	$job_info_file = file(get_path("sample_status_folder")."/jobinfo.txt");
	$ps_name_to_jobIDs=array();
	foreach($job_info_file as $line)//build a row for each entry in jobinfo.txt
	{
		list($jobnumber,,,,$found_proc_sample_name)=explode("\t",$line);
		if ($found_proc_sample_name==$trio)
		{
			$ps_name_to_jobIDs[]=$jobnumber;
		}
	}
	//check if a job corresponding to the family is still in queue
	foreach($ps_name_to_jobIDs as $jobID)
	{
		$qstat_output="";
		$qstat_output=shell_exec("qstat -j ".$jobID." 2> /dev/null"); //don't show error messages
		if ($qstat_output!="")//if the job was found
		{
			trigger_error("Family ".$trio. " is already in queue with jobnumber=".$jobID."!", E_USER_ERROR);
		}
	}
}

//find the project name and process sample type for a sample ID, find and check runnumber
function get_parameters($processed_sample_name)
{

	//get sample name and process ID
	list($sample_name,$process_id) = explode("_",$processed_sample_name);
	$process_id=ltrim($process_id,'0');
	$db_connect = DB::getInstance('NGSD');
	
	//get sample id
	$result = $db_connect->executeQuery("SELECT id FROM  sample WHERE name='".$sample_name."'");
	$sample_id=$result[0]["id"];
	
	//get processed processed sample id and project id
	$result = $db_connect->executeQuery("SELECT id,project_id FROM  processed_sample WHERE sample_id=".$sample_id." AND process_id= ".$process_id);
	$proc_sample_id=$result[0]["id"];
	$projectID=$result[0]["project_id"];
	
	//get project name
	$result = $db_connect->executeQuery("SELECT name,type FROM project WHERE id=".$projectID);
	$projectname=$result[0]["name"];
	$project_type=$result[0]["type"]; 

	//determine project folder
	$project_folder = get_path("project_folder")."/".$project_type."/".$projectname."/";
	$sample_bam = $project_folder."/Sample_".$processed_sample_name."/".$processed_sample_name.".bam";
	if(!file_exists($sample_bam))	trigger_error("Could not find sample bam '$sample_bam'", E_USER_ERROR);
	
	return array("ps_id" => $proc_sample_id, "project_id" => $projectID, "process_id"=>$process_id, "project_type"=>$project_type, "project_name"=>$projectname,"path"=>$sample_bam, "wd"=>$project_folder);
}

check_trio_in_queue($f,$m,$c);

$f_params=get_parameters($f);
$m_params=get_parameters($m);
$c_params=get_parameters($c);

//create output folder
$out_folder = $c_params["wd"]."/Trio_".$c."_".$f."_".$m."/";
if (!file_exists($out_folder)) mkdir($out_folder);

//queue trio for analysis
$command = "php ".repository_basedir()."/src/Pipelines/trio.php";
$args = "-c ".$c_params["path"]." -f ".$f_params["path"]." -m ".$m_params["path"]." -out_folder ".$out_folder." --log ".$out_folder."trio.log";
$sample_status = get_path("sample_status_folder")."/data/";
$queue = $high_priority ? "-q NGSlong,re_analysis,srv016_long" : "-q NGSlong,srv016_long";
$queue_command = "qsub -V -b y -wd /tmp/ -m n -M florian.lenz@med.uni-tuebingen.de -e $sample_status $queue -o $sample_status $command $args";
$qsub_return_line=shell_exec($queue_command);

//extract job number
$exploded_qsub_return_line = explode(" ", $qsub_return_line);
$jobnumber = intval($exploded_qsub_return_line[2]);
if($jobnumber>=1)
{
	$user_name = trim(shell_exec('whoami'));
	$outputline = array($jobnumber,date("d-m-Y_H:i:s"), "n/a", "n/a", $c."|".$f."|".$m, $c_params["ps_id"]."|".$f_params["ps_id"]."|".$m_params["ps_id"], "n/a", "n/a", $user_name);
	$outputline[] = $high_priority ? "high" : "low";
	file_put_contents(get_path("sample_status_folder")."/jobinfo.txt", implode("\t", $outputline)."\n", FILE_APPEND);
}
else
{
	trigger_error("Could not queue job - job number lower than 1 was returned!", E_USER_ERROR);
}