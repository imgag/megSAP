<?php 
/** 
	@page queue_sample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");


$parser = new ToolBase("queue_sample", "Queues sample for analysis.");
$parser->addString("sample",  "Processed sample identifier.", false);
$parser->addString("steps", "Comma-separated list of processing steps to perform.", true, "ma,vc,an,db,cn");
//optional
$parser->addString("user", "Name of the user who queued the analysis (current user if unset).", true, "");
$parser->addFlag("high_priority", "Assign a high priority to job.");
extract($parser->parse($argv));

//init
if ($user=="")
{
	$user = exec('whoami');
}

//check if the given processed sample is already in queue
$jobinfo = sge_jobinfo();
foreach($jobinfo as $id => $data)
{
	$command = $data[1];
	if (!contains($command, "analyze.php")) continue;
	if (contains($command, "-name {$sample}"))	
	{
		trigger_error("Analysis already running with ID $id", E_USER_ERROR);
	}
}

//get information from NGSD 
$info = get_processed_sample_info($sample);
$project_folder = $info["project_folder"];

//check sample folder is available
$sample_folder = $project_folder."Sample_".$sample."/";
if(!is_dir($sample_folder))
{
	trigger_error("Could not find sample folder '$sample_folder'!", E_USER_ERROR);
}

//check that fastq files are there
$files = glob($sample_folder.$sample."*.fastq.gz");
if(count($files)<2 && $info['sys_type'] != "RNA")
{
	trigger_error("Could not find at least two FASTQ files starting with '$sample' in '$sample_folder'!", E_USER_ERROR);
}

//determine command and arguments
if($info['is_tumor'] && $info['normal_name']!="" && $info['sys_type'] != "RNA")
{	
	$outfolder = $project_folder."/Somatic_".$sample."-".$info['normal_name']."/";
	if (!file_exists($outfolder))
	{
		mkdir($outfolder);		
		if (!chmod($outfolder, 0777))
		{
			trigger_error("Could not change privileges of folder '{$outfolder}'!", E_USER_ERROR);
		}
	}
	
	//determine somatic steps
	$steps_som = array_intersect(array("ma", "vc", "an", "ci", "db"), explode(",",$steps));
	if (in_array("an", $steps_som))
	{
		$steps_som[] = "ci"; 
	}
	$steps_som = implode(",",$steps_som);
	
	$command = "php ".repository_basedir()."/src/Pipelines/somatic_dna.php";
	$args = "-p_folder {$project_folder} -t_id {$sample} -n_id ".$info['normal_name']." -o_folder {$outfolder} -steps {$steps_som} --log {$outfolder}somatic_dna_".date("YmdHis").".log";
}
elseif ($info['sys_type'] == "RNA")
{
	$command = "php ".repository_basedir()."/src/Pipelines/analyze_rna.php";
	
	//if steps argument is default, replace with analyze_rna default value
	if ($steps == "ma,vc,an,db,cn")
	{
		$steps = "ma,rc,an,fu,db";
	}
	else
	{
		//reduce to valid steps for analyze_rna
		$steps = implode(",", array_intersect(explode(",", $steps), explode(",", "ma,rc,an,fu,db")));
	}
	$args = "-folder {$sample_folder} -name {$sample} -steps {$steps} --log {$sample_folder}analyze_rna_".date("YmdHis").".log";
}
else
{
	$command = "php ".repository_basedir()."/src/Pipelines/analyze.php";
	$args = "-folder {$sample_folder} -name {$sample} -steps {$steps} --log {$sample_folder}analyze_".date("YmdHis").".log";
	if ($info['sys_type']=="WGS") //whole genome => use 5 threads (default is 2)
	{
		$args .= " -threads 5";
	}	
}

//submit to queue
$queues = explode(",", get_path("queues_default"));
if($high_priority)
{
	$queues = array_merge($queues, explode(",", get_path("queues_high_priority")));
}
elseif ($info['sys_type']=="RNA")
{
	$queues = explode(",", get_path("queues_high_mem"));
}

$sample_status = get_path("sample_status_folder")."/data/";
$slots = $info['sys_type']=="WGS" ? "-pe smp 2" : ""; //use two instead of one slots for WGS

$qsub_return_line = exec("qsub -V $slots -b y -wd $project_folder -m n -M ".get_path("queue_email")." -e $sample_status -q ".implode(",", $queues)." -o $sample_status $command $args");

//extract job number
$exploded_qsub_return_line = explode(" ", $qsub_return_line);
$jobnumber = intval($exploded_qsub_return_line[2]);
if($jobnumber>0)
{
	$outputline = array($jobnumber, date("d-m-Y_H:i:s"), $info['run_name'], $info['run_id'], $sample, $info['ps_id'], $info['project_name'], $info['project_id'], $user, $high_priority ? "high" : "low");
	file_put_contents(get_path("sample_status_folder")."/jobinfo.txt", implode("\t", $outputline)."\n", FILE_APPEND);
}
else
{
	trigger_error("Could not queue job - job number lower than 1 was returned!", E_USER_ERROR);
}