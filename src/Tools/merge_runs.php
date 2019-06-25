<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$parser = new ToolBase("merge_runs", "Creates Makefile from two runs for merging single samples.");
$parser->addString("run","Run-ID of run containing samples where other samples are merged into.",false);
$parser->addString("run_old","ID of run that contains samples which are merged into samples of other run.",false);
$parser->addString("merge_script","path to script which is used to merge samples.",true,realpath(dirname($_SERVER['SCRIPT_FILENAME'])."/../Tools/merge_samples.php"));
$parser->addString("out","Path to output Makefile. Writes to stdout otherwise.",true);
extract($parser->parse($argv));

$db = DB::getInstance("NGSD");

//parse sample ids from run
$query = "SELECT run.id, ps.sample_id, ps.process_id, s.name FROM sequencing_run as run, processed_sample as ps, sample as s WHERE run.name = '{$run}' AND ps.sequencing_run_id = run.id AND s.id = ps.sample_id";
$res = $db->executeQuery($query);
$sample_names = array();
$sample_names_old = array();
foreach($res as $data)
{
	//Include only samples that are included in the old run to be merged
	$processed_id_old = $data["process_id"]-1;
	
	$matching_old_sample_found = false;
	
	//Take "youngest" precursor
	while($processed_id_old > 0)
	{
		$temp_sample_name_old = $data["name"] . "_" . str_pad($processed_id_old,2,'0',STR_PAD_LEFT);
		
		$info_old = get_processed_sample_info($db,$temp_sample_name_old);
		if($info_old["run_name"] == $run_old)
		{
			$matching_old_sample_found = true;
			break;
		}
		--$processed_id_old;
	}
	
	if($matching_old_sample_found)
	{
		$sample_names[] = $data["name"] . "_" . str_pad($data["process_id"],2,'0',STR_PAD_LEFT);
		$sample_names_old[] = $data["name"] . "_" . str_pad($processed_id_old,2,'0',STR_PAD_LEFT);
	}
}

$out_text = "all: help merge\n";
$out_text .= "help:\n";
$out_text .= "\t@echo Makefile was created using " . __FILE__ . "\n";
$out_text .= "\t@echo All merge commands should be checked manually before merging.\n";
$out_text .= "merge:\n";

for($i=0;$i<count($sample_names);++$i)
{	
	$info_new = get_processed_sample_info($db,$sample_names[$i]);
	$info_old = get_processed_sample_info($db,$sample_names_old[$i]);
	
	//Skip samples with differing target regions
	if($info_new["sys_target"] != $info_old["sys_target"])
	{
		trigger_error("Target region  of ". $sample_names[$i] . " and " . $sample_names_old[$i] . " is not the same. Aborting these samples.\n",E_USER_WARNING);
		continue;
	}

	//Skip samples not included in old run
	if(isset($run_old))
	{
		if($info_old["run_name"] != $run_old) continue;
	}
	
	$out_text .= "\tphp {$merge_script} -into " .$sample_names[$i] . " -ps ". $sample_names_old[$i] ."\n";
}
if(isset($out))
{
	file_put_contents($out,$out_text);
}
else
{
	print($out_text);
}
?>