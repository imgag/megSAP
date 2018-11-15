<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("reannotate_unsolved_diagnostic_cases", "Re-annotation of unsolved diagnostic cases (not bad quality, WES/WGS).");
$parser->addInt("sequenced_ago", "Maximum number of days the sample was sequenced ago.", true, 365); //1 year
$parser->addInt("annotated_ago", "Maximum number of days the sample was last annotated ago.", true, 42); //6 weeks
extract($parser->parse($argv));

//select diagnostic samples with good/medium quality that were sequenced in the defined time range
$db = DB::getInstance("NGSD");
$query = "SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) as id ".
			 "FROM processed_sample ps, sample s, project p, sequencing_run r, diag_status ds, processing_system sys ".
			 "WHERE sys.id=ps.processing_system_id AND ps.sample_id = s.id AND ps.project_id=p.id AND ds.processed_sample_id=ps.id AND (sys.type='WGS' OR sys.type='WES') AND p.type='diagnostic' AND ps.sequencing_run_id=r.id AND ps.quality!='bad' AND ds.outcome='n/a' AND r.start_date>SUBDATE(NOW(), INTERVAL $sequenced_ago DAY)";
$result = $db->executeQuery($query);

//construct date
$anno_date = new DateTime();
$anno_date->sub(new DateInterval("P{$annotated_ago}D"));
$anno_timestamp = $anno_date->getTimestamp();

//re-queue annotation
for ($i=0; $i<count($result); ++$i)
{
	$ps = $result[$i]["id"];
	$info = get_processed_sample_info($db, $ps, true);
	print ($i+1)."/".count($result)." {$ps} (".$info['project_name'].") ";
	
	//skip when folder missing
	if (!file_exists($info['ps_folder']))
	{
		print "Skipped: Sample folder missing!\n";
		continue;
	}
	
	//skip when BAM missing
	$gsvar = substr($info['ps_bam'],0, -4).".GSvar";
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
	
	//queue
	print "Queuing...\n";
	list($stdout, $stderr, $return) = $parser->execTool("NGS/db_queue_analysis.php", "-type 'single sample' -samples $ps -args '-steps an'", false);
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