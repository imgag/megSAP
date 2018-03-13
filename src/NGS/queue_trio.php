<?php 
/** 
	@page queue_trio
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("queue_trio", "Queues a trio for analysis.");
$parser->addString("c",  "Processed sample identifier of child", false);
$parser->addString("f",  "Processed sample identifier of father", false);
$parser->addString("m",  "Processed sample identifier of mother", false);
//optional
$parser->addString("start", "Start step.", true, "check");
$parser->addString("user", "Name of the user who queued the analysis (current user if unset).", true, "");
$parser->addFlag("high_priority", "Assign a high priority to job.");
extract($parser->parse($argv));

//init
if ($user=="")
{
	$user = exec('whoami');
}

//check if a trio analysis for that family is already in queue
$jobinfo = sge_jobinfo();
foreach($jobinfo as $id => $data)
{
	$command = $data[1];
	if (!contains($command, "trio.php")) continue;
	if (contains($command, $f) && contains($command, $m) && contains($command, $c))
	{
		trigger_error("Analysis already running with ID $id", E_USER_ERROR);
	}
}

$f_params = get_processed_sample_info($f);
$m_params = get_processed_sample_info($m);
$c_params = get_processed_sample_info($c);

//create output folder
$out_folder = $c_params["project_folder"]."/Trio_".$c."_".$f."_".$m."/";
if (!file_exists($out_folder))
{
	mkdir($out_folder);
	if (!chmod($out_folder, 0777))
	{
		trigger_error("Could not change privileges of folder '{$out_folder}'!", E_USER_ERROR);
	}
}

//determine command and arguments
$command = "php ".repository_basedir()."/src/Pipelines/trio.php";
$args = "-c ".$c_params["ps_bam"]." -f ".$f_params["ps_bam"]." -m ".$m_params["ps_bam"]." -out_folder {$out_folder} --log {$out_folder}trio.log";

//queue trio for analysis
$queues = explode(",", get_path("queues_default"));
if($high_priority)
{
	$queues = array_merge($queues, explode(",", get_path("queues_high_priority")));
}
$sample_status = get_path("sample_status_folder")."/data/";
$qsub_return_line = exec("qsub -V -b y -wd /tmp/ -m n -M ".get_path("queue_email")." -e $sample_status -q ".implode(",", $queues)." -o $sample_status $command $args");

//extract job number
$exploded_qsub_return_line = explode(" ", $qsub_return_line);
$jobnumber = intval($exploded_qsub_return_line[2]);
if($jobnumber>=1)
{
	$outputline = array($jobnumber, date("d-m-Y_H:i:s"), "n/a", "n/a", "trio|".$c."|".$f."|".$m, $c_params["ps_id"]."|".$f_params["ps_id"]."|".$m_params["ps_id"], "n/a", "n/a", $user, $high_priority ? "high" : "low");
	file_put_contents(get_path("sample_status_folder")."/jobinfo.txt", implode("\t", $outputline)."\n", FILE_APPEND);
}
else
{
	trigger_error("Could not queue job - job number lower than 1 was returned!", E_USER_ERROR);
}