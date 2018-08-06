<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("batch_reanalysis", "Batch re-analysis of single-sample analyses.");
$parser->addStringArray("samples", "Processed sample names.", false);
$parser->addString("steps", "Analysis steps to perform.", true);
$parser->addInt("threads", "The maximum number of threads used.", true);
$parser->addString("before", "Only samples analyzed before the given date in format DD.MM.YYYY are reanalyzed.", true);
extract($parser->parse($argv));

//convert 'before' to timestamp
if (isset($before))
{
	$converted = DateTime::createFromFormat("d.m.Y", $before);
	if(!$converted)
	{
		trigger_error("Could not convert '$before' to a date!", E_USER_ERROR);
	}
	$before = $converted->getTimestamp();
}

$count = count($samples);
for ($i=1; $i<=$count; ++$i)
{
	$ps = trim($samples[$i-1]);
	
	//determine folder name
	$db_conn = DB::getInstance("NGSD");
	$info = get_processed_sample_info($db_conn, $ps, false);
	if (is_null($info))
	{
		print "$i/$count: Skipped '$ps' because it was not found in NGSD!\n";
		continue;
	}
	$folder = $info['ps_folder'];
	
	//check folder exists
	if(!file_exists($folder))
	{
		print "$i/$count: Skipped '$ps' because folder does not exist: $folder\n";
		continue;
	}
	
	//check last analysis date 
	if (isset($before))
	{
		$bam = $info['ps_bam'];
		$gsvar = substr($bam, 0, -3)."GSvar";
		if (file_exists($bam) && filemtime($bam)>$before && file_exists($gsvar) && filemtime($gsvar)>$before)
		{
			print "$i/$count: Skipped '$ps' because it was analyzed after ".date("m.d.Y", $before)."\n";
			continue;
		}
	}
	
	//variant calling with multi-sample pipeline
	print "$i/$count: Analyzing '$folder'...\n";
	$args = array();
	$args[] = "-name $ps";
	$args[] = "-folder $folder";
	if (isset($steps))
	{
		$args[] = "-steps $steps";
	}
	if (isset($threads))
	{
		$args[] = "-threads $threads";
	}
	list($stdout, $stderr, $return) = $parser->execTool("Pipelines/analyze.php", implode(" ", $args), false);
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