<?php 
/** 
	@page queue_multi
	
	@todo check if analysis is already in queue
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("queue_multi", "Adds a multi-sample analysis to the queue.");
$parser->addStringArray("samples", "Processed sample names files.", false);
$parser->addStringArray("status", "Affected status of the samples - 'affected' or 'control'.", false);
$parser->addFlag("high_priority", "Assign a high priority to the job.");
extract($parser->parse($argv));

//get BAM files for samples
$bams = array();
$db = DB::getInstance('NGSD');
foreach($samples as $sample)
{
	list($sample_name, $process_id) = explode("_", $sample);
	$res = $db->executeQuery("SELECT p.name, p.type FROM project p, processed_sample ps, sample s WHERE ps.project_id=p.id AND ps.sample_id=s.id AND s.name='$sample_name' AND ps.process_id='".(int)$process_id."'");
	if (count($res)!=1)
	{
		trigger_error("Could not determine project information for sample '$sample'!", E_USER_ERROR);
	}
	$bams[] = get_path("project_folder")."/".$res[0]['type']."/".$res[0]['name']."/Sample_".$sample."/".$sample.".bam";
}

//create output folder (in project folder of first sample)
$out_folder = dirname(dirname($bams[0]))."/Multi_".implode("_", $samples);
if (!file_exists($out_folder)) 
{
	mkdir($out_folder);
}

//determine command and arguments
$command = "php ".repository_basedir()."/src/Pipelines/multisample.php";
$args = "-bams ".implode(" ", $bams)." -status ".implode(" ", $status)." -out_folder ".$out_folder." --log ".$out_folder."multi.log";

//queue trio for analysis
$queues = explode(",", get_path("queues_default"));
if($high_priority)
{
	$queues = array_merge($queues, explode(",", get_path("queues_high_priority")));
}
$sample_status = get_path("sample_status_folder")."/data/";
exec2("qsub -V -b y -wd /tmp/ -m n -M ".get_path("queue_email")." -e $sample_status -q ".implode(",", $queues)." -o $sample_status $command $args");
