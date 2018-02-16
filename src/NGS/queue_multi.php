<?php 
/** 
	@page queue_multi
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("queue_multi", "Adds a multi-sample analysis to the queue.");
$parser->addStringArray("samples", "Processed sample names files.", false);
$parser->addStringArray("status", "Affected status of the samples - 'affected' or 'control'.", false);
//optional
$parser->addString("steps", "Processing steps to perform.", true);
$parser->addString("user", "Name of the user who queued the analysis (current user if unset).", true, "");
$parser->addFlag("high_priority", "Assign a high priority to the job.");
extract($parser->parse($argv));

//init
if ($user=="")
{
	$user = exec('whoami');
}

//check if analysis is already running
$jobinfo = sge_jobinfo();
foreach($jobinfo as $id => $data)
{
	$command = $data[1];
	if (!contains($command, "multisample.php")) continue;
	$all_found = true;
	foreach($samples as $sample)
	{
		if (!contains($command, $sample)) $all_found = false;
	}
	if ($all_found)
	{
		trigger_error("Analysis already running with ID $id", E_USER_ERROR);
	}
}

//get information from NGSD
$ps_ids = array();
$bams = array();
$project_folder = "";
foreach($samples as $sample)
{
	$info = get_processed_sample_info($sample);
	$ps_ids[] = $info['ps_id'];
	$bams[] = $info['ps_bam'];
	if ($project_folder=="")
	{
		$project_folder = $info['project_folder'];
	}
}

//create output folder
sort($samples);
$out_folder = $project_folder."/Multi_".implode("_", $samples)."/";
if (!file_exists($out_folder)) 
{
	mkdir($out_folder);
	if (!chmod($out_folder, 0777))
	{
		trigger_error("Could not change privileges of folder '{$out_folder}'!", E_USER_ERROR);
	}
}

//determine command and arguments
$command = "php ".repository_basedir()."/src/Pipelines/multisample.php";
$args = array("-bams ".implode(" ", $bams), "-status ".implode(" ", $status), "-out_folder {$out_folder}", "--log {$out_folder}multi.log");
if ($steps!="")
{
	$args[]  = "-steps {$steps}"; 
}

//queue trio for analysis
$queues = explode(",", get_path("queues_default"));
if($high_priority)
{
	$queues = array_merge($queues, explode(",", get_path("queues_high_priority")));
}
$sample_status = get_path("sample_status_folder")."/data/";
$qsub_return_line = exec("qsub -V -b y -wd /tmp/ -m n -M ".get_path("queue_email")." -e $sample_status -q ".implode(",", $queues)." -o $sample_status $command ".implode(" ", $args));

//extract job number
$exploded_qsub_return_line = explode(" ", $qsub_return_line);
$jobnumber = intval($exploded_qsub_return_line[2]);
if($jobnumber>=1)
{
	$outputline = array($jobnumber, date("d-m-Y_H:i:s"), "n/a", "n/a", "multi|".implode("|", $samples)."|".implode("|", $status), implode("|", $ps_ids), "n/a", "n/a", $user, $high_priority ? "high" : "low");
	file_put_contents(get_path("sample_status_folder")."/jobinfo.txt", implode("\t", $outputline)."\n", FILE_APPEND);
}
else
{
	trigger_error("Could not queue job - job number lower than 1 was returned!", E_USER_ERROR);
}
