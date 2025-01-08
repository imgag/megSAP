<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("batch_reanalysis_multi", "Batch re-analysis of multi-sample analyses.");
$parser->addStringArray("files", "Mutli-sample GSvar files.", false);
$parser->addString("steps", "Analysis steps to perform.", false);
$parser->addInt("threads", "Number of threads used.", true);
$parser->addString("before", "Only samples analyzed before the date are reanalyzed (checks GSvar file date, format 'DD.MM.YYYY').", true);
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

$count = count($files);
for ($i=1; $i<=$count; ++$i)
{
	$gsvar = trim($files[$i-1]);
	
	//check last analysis date 
	if (filemtime($gsvar)>$before)
	{
		print "$i/$count: Skipped '$gsvar' because it was analyzed after ".date("m.d.Y", $before)."\n";
		continue;
	}
	
	//extract sample status from headers
	$ps2status = array();
	list($headers) = exec2("grep '##SAMPLE=' $gsvar", false);
	foreach($headers as $header)
	{
		$ps = null;
		$parts = explode(",", substr(trim($header), 10, -1));
		foreach($parts as $p)
		{
			@list($name, $value) = explode("=", $p);
			if ($name=="ID") $ps = $value;
			if ($name=="DiseaseStatus") $ps2status[$ps] = $value;
		}
	}
	
	//determine BAM/status 
	$bams = array();
	$status = array();
	$out_folder = dirname($gsvar);
	$project_folder = dirname($out_folder);
	$tmp = explode("_", basename($out_folder));
	for($j=1; $j<count($tmp); $j+=2)
	{
		$ps = $tmp[$j]."_".$tmp[$j+1];
		$bam = $project_folder."/Sample_{$ps}/{$ps}.bam";
		if (!file_exists($bam))
		{
			print "$i/$count: Skipped '$out_folder' because BAM file not found at: $bam\n";
			continue 2;
		}
		$bams[] = $bam;
		if (!isset($ps2status[$ps]))
		{
			print "$i/$count: Skipped '$out_folder' because sample status for $ps could not be determined from header of $gsvar\n";
			continue 2;
		}
		$status[] = $ps2status[$ps];
	}
	
	//perform analysis
	print "$i/$count: Analyzing '$gsvar'...\n";
	
	$args = array();
	$args[] = "-bams ".implode(" ", $bams);
	$args[] = "-status ".implode(" ", $status);
	$args[] = "-out_folder {$out_folder}";
	$args[] = "-steps {$steps}";
	$args[] = "--log {$out_folder}/multi.log";
	if (isset($threads))
	{
		$args[] = "-threads {$threads}";
	}
	list($stdout, $stderr, $return) = $parser->execTool("Pipelines/multisample.php", implode(" ", $args), false);
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