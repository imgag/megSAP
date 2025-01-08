<?php 
/** 
	@page db_queue_analysis
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("db_queue_analysis", "Queues an analysis in the NGSD.");
$parser->addEnum("type",  "Analysis type.", false, array('single sample','multi sample','trio','somatic'));
$parser->addStringArray("samples", "Processed sample name(s).", false);
//optional
$parser->addStringArray("info", "Sample info entries for complex analysis jobs: 'affected/control' for multi sample, 'tumor/normal/tumor_rna' for somatic, 'child/father/mother' for trio.", true);
$parser->addFlag("high_priority", "Perform analysis with high priority");
$parser->addString("args",  "Custom arguments passed on to the analysis script.", true);
$parser->addString("user", "Name of the user who queued the analysis (current user if unset).", true, "");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//check info
if ($type=="single sample") //single sample => no info
{
	$info = array_fill(0, count($samples), "");
}
if (count($samples)!=count($info))
{
	trigger_error("Mismatching counts of samples (".count($samples).") and status (".count($info).")!", E_USER_ERROR);
}

//check args
if (is_null($args)) $args = "";

//get user ID from NGSD
if ($user=="") $user = exec('whoami');
$db_conn = DB::getInstance($db);
$user_id = $db_conn->getValue("SELECT id FROM user WHERE user_id='".$user."' AND active='1'", -1);
if ($user_id==-1)
{
	trigger_error("User '$user' not found in NGSD!", E_USER_ERROR);
}

//get processed sample data from NGSD
$ps_ids = array();
foreach($samples as $sample)
{
	$ps_ids[] = get_processed_sample_id($db_conn, $sample);
}

//check if the analysis is already in queue
$result = $db_conn->executeQuery("SELECT j.id FROM analysis_job j, analysis_job_sample js WHERE j.type=:type AND js.analysis_job_id=j.id AND js.processed_sample_id=:ps_id ORDER BY id DESC", array("type"=>$type, "ps_id"=>$ps_ids[0]));
foreach($result as $row)
{
	$job_id = $row['id'];
	$job_info = analysis_job_info($db_conn, $job_id);
		
	//only running jobs
	$last_status = end($job_info['history']);
	if ($last_status!="queued" && $last_status!="started") continue;
	
	//same samples
	if (count($job_info['samples'])!=count($samples)) continue;
	$same_samples = true;
	for ($i=0; $i<count($samples); ++$i)
	{
		if ($job_info['samples'][$i]!=$samples[$i]."/".$info[$i])
		{
			$same_samples = false;
			break;
		}
	}
	if (!$same_samples) continue;
	
	trigger_error("Analysis is already running with NGSD job ID '{$job_id}' and last status '{$last_status}'!", E_USER_ERROR);
}

//queue analysis
$db_conn->beginTransaction();
$db_conn->executeStmt("INSERT INTO `analysis_job`(`type`, `high_priority`, `args`) VALUES (:type, :prio, :args)", array("type"=>$type, "prio"=>($high_priority ? "1" : "0"), "args"=>$args));
$job_id = $db_conn->lastInsertId();
$db_conn->executeStmt("INSERT INTO `analysis_job_history`(`analysis_job_id`, `time`, `user_id`, `status`, `output`) VALUES ({$job_id}, '".get_timestamp(false)."', $user_id, 'queued', '')");
for($i=0; $i<count($ps_ids); ++$i)
{
	$db_conn->executeStmt("INSERT INTO `analysis_job_sample`(`analysis_job_id`, `processed_sample_id`, `info`) VALUES ({$job_id},{$ps_ids[$i]},:info)", array("info"=>$info[$i]));
}
$db_conn->endTransaction();

?>
