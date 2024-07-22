<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("reannotate_unsolved_diagnostic_cases", "Re-annotation of unsolved diagnostic cases (not bad quality, WES/WGS).");
$parser->addInt("sequenced_ago", "Maximum number of days the sample was sequenced ago.", true, 240); //~ 8 months
$parser->addInt("annotated_ago", "Maximum number of days the sample was last annotated ago.", true, 35); //5 weeks
extract($parser->parse($argv));

//select diagnostic samples with good/medium quality that were sequenced in the defined time range
$db = DB::getInstance("NGSD");
$query = "SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) as id, ds.outcome, sys.type ".
			 "FROM sample s, project p, sequencing_run r, processing_system sys, processed_sample ps LEFT JOIN diag_status ds ON ds.processed_sample_id=ps.id ".
			 "WHERE sys.id=ps.processing_system_id AND ps.sample_id=s.id AND s.disease_status!='Unaffected' AND ps.project_id=p.id AND (sys.type='WGS' OR sys.type='WES' OR sys.type='lrGS') AND p.type='diagnostic' AND s.tumor='0' AND ps.sequencing_run_id=r.id AND ps.quality!='bad' AND r.start_date>SUBDATE(NOW(), INTERVAL $sequenced_ago DAY) AND r.end_date<NOW() AND ps.id NOT IN (SELECT processed_sample_id FROM merged_processed_samples)";
$result = $db->executeQuery($query);

//skip if outcome is set
$tmp = array();
for ($i=0; $i<count($result); ++$i)
{
	$outcome = $result[$i]["outcome"];
	if ($outcome!="n/a" && $outcome!="") continue;
	
	$tmp[] = $result[$i];
}
$result = $tmp;

//construct date
$anno_date = new DateTime();
$anno_date->sub(new DateInterval("P{$annotated_ago}D"));
$anno_timestamp = $anno_date->getTimestamp();

//re-queue annotation
for ($i=0; $i<count($result); ++$i)
{
	$ps = $result[$i]["id"];
	$info = get_processed_sample_info($db, $ps);
	print ($i+1)."/".count($result)." {$ps} (".$info['project_name'].") ";
	
	//skip when folder missing
	if (!file_exists($info['ps_folder']))
	{
		print "Skipped: Sample folder missing!\n";
		continue;
	}
	
	//skip when GSvar missing
	$gsvar = $info['ps_folder']."/".$ps.".GSvar";
	if (!file_exists($gsvar))
	{
		print "Skipped: GSvar file missing!\n";
		continue;
	}
	
	//skip when too new
	if (filemtime($gsvar)>=$anno_timestamp)
	{	
		print "Skipped: was analyzed after ".$anno_date->format("d.m.Y")."\n";
		continue;
	}
	
	//define pipeline arguments based on system type
	$sys_type = $result[$i]["type"];
	$args = ($sys_type=="lrGS") ? "-steps an" : "-steps vc,cn,sv -annotation_only";
	
	//queue
	print "Queuing...\n";
	list($stdout, $stderr, $return) = $parser->execTool("NGS/db_queue_analysis.php", "-type 'single sample' -samples $ps -args '{$args}'", false);
	if ($return!=0)
	{
		print "  Error occurred:\n";
		$lines = array_merge($stdout, $stderr);
		foreach($lines as $line)
		{
			print "  ".trim($line)."\n";
		}
	}
}


?>
