<?php
/** 
	@page find_wrong_adapters
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//extracts adapter sequence from a log file line
function extract_adapter($input)
{
	$parts = explode(":", $input);
	return trim(end($parts));
}

//determines the percentage of mismatches between the first n bases of two sequences (N is counted as match)
function mm_perc($a, $c, $n)
{
	$mm = 0;
	for ($i=0; $i<$n; ++$i)
	{
		if ($c[$i]=="N") continue;
		$mm += $a[$i]!=$c[$i];
	}
	return $mm/$n;
}

//parse command line arguments
$parser = new ToolBase("find_wrong_adapters", "Checks for mismatches between given and consensus adapter in all samples.");
extract($parser->parse($argv));

//get list of mapping qcML files
$project_folders  = get_path("project_folder");
list($qc_files) = exec2("find -L ".$project_folders['diagnostic']." ".$project_folders['research']."  -maxdepth 3 -name \"*_stats_map.qcML\"");

//NGSD connection
$db = DB::getInstance("NGSD");

print "#sample\tsystem\ta1_len\ta1_mismatch_perc\ta2_len\ta2_mismatch_perc\n";
foreach($qc_files as $file)
{
	//check sample exists in NGSD
	list($name, $ps_id) = explode("_", str_replace("_stats_map.qcML", "", basename($file))."_");
	$s_id = $db->getID("sample", "name", $name, false);
	if($s_id==-1) continue;
	
	//check processed sample exists in NGSD and extract processing system name
	$res = $db->executeQuery("SELECT ps.id, sys.name_manufacturer FROM processed_sample ps, processing_system sys WHERE ps.processing_system_id=sys.id AND ps.sample_id='$s_id' AND ps.process_id='".(int)($ps_id)."'");
	if(count($res)!=1) continue;
	$sys = $res[0]['name_manufacturer'];
	
	//grep for SeqPurge output
	list($matches, $stderr) = exec2("grep 'adapter sequence (' $file", false);
	if(trim(implode("", $stderr))!="") trigger_error("Grep error: ".implode("\n", $stderr), E_USER_ERROR);
	
	//skip log files that contain no/several SeqPurge outputs
	if(count($matches)!=4) continue;
	
	//extract adapters
	list($a1, $a1_con, $a2, $a2_con) = array_map("extract_adapter", $matches);	
	
	//calculate length/distance
	$a1_len = min(strlen($a1), strlen($a1_con));
	if ($a1_len<10) continue;
	$a1_mmp = mm_perc($a1, $a1_con, $a1_len);

	//calculate length/distance
	$a2_len = min(strlen($a2), strlen($a2_con));
	if ($a2_len<10) continue;
	$a2_mmp = mm_perc($a2, $a2_con, $a2_len);
	
	print "{$name}_{$ps_id}\t$sys\t$a1_len\t$a1_mmp\t$a2_len\t$a2_mmp\n";
}

?>