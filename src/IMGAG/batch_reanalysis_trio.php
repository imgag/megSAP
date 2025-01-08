<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("batch_reanalysis_trio", "Batch re-analysis of trio analyses.");
$parser->addStringArray("files", "Trio GSvar files.", false);
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
	
	//determine sample BAM files
	$out_folder = dirname($gsvar);
	$project_folder = dirname($out_folder);
	$tmp = explode("_", basename($out_folder));
	$ps_c = $tmp[1]."_".$tmp[2];
	$ps_f = $tmp[3]."_".$tmp[4];
	$ps_m = $tmp[5]."_".$tmp[6];
	$bam_c = $project_folder."/Sample_{$ps_c}/{$ps_c}.bam";	
	if(!file_exists($bam_c))
	{
		print "$i/$count: Skipped '$out_folder' because child BAM file not found at: $bam_c\n";
		continue;
	}
	$bam_f = $project_folder."/Sample_{$ps_f}/{$ps_f}.bam";
	if(!file_exists($bam_f))
	{
		print "$i/$count: Skipped '$out_folder' because father BAM file not found at: $bam_f\n";
		continue;
	}
	$bam_m = $project_folder."/Sample_{$ps_m}/{$ps_m}.bam";
	if(!file_exists($bam_m))
	{
		print "$i/$count: Skipped '$out_folder' because mother BAM file not found at: $bam_m\n";
		continue;
	}
	
	//perform analysis
	print "$i/$count: Analyzing '$gsvar'...\n";
	
	$args = array();
	$args[] = "-f {$bam_c}";
	$args[] = "-m {$bam_f}";
	$args[] = "-c {$bam_m}";
	$args[] = "-out_folder {$out_folder}";
	$args[] = "-steps {$steps}";
	$args[] = "--log {$out_folder}/trio.log";
	if (isset($threads))
	{
		$args[] = "-threads {$threads}";
	}
	list($stdout, $stderr, $return) = $parser->execTool("Pipelines/trio.php", implode(" ", $args), false);
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