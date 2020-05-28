<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("batch_reanalysis", "Batch re-analysis of single-sample analyses.");
$parser->addStringArray("samples", "Processed sample names.", false);
$parser->addString("steps", "Analysis steps to perform.", false);
$parser->addEnum("mode", "Excution mode: 'default' executes the analysis sequentially in this script, 'print' only prints samples, but performs no analysis, 'sge' queues the analysis in SGE.", true, array("default", "print", "sge"), "default");
$parser->addString("before", "Only samples analyzed before the date are reanalyzed (considers 'steps', format 'DD.MM.YYYY').", true);
$parser->addInt("threads", "Number of threads used (for 'default' mode).", true);
$parser->addInt("max", "Maximum number of jobs to start (for 'sge' mode).", true);
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

$analyzed = 0;
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
		$skip = true;
		
		//check BAM
		$bam = $info['ps_bam'];
		if (contains($steps, "ma"))
		{
			if (!file_exists($bam) || filemtime($bam)<$before)
			{
				$skip = false;
			}
		}
		
		//check VCF
		if (contains($steps, "vc"))
		{
			$vcf = substr($bam, 0, -4)."_var.vcf.gz";
			if (!file_exists($vcf) || filemtime($vcf)<$before)
			{
				$skip = false;
			}
		}
		
		//check GSvar
		if (contains($steps, "an"))
		{
			$gsvar = substr($bam, 0, -3)."GSvar";
			if (!file_exists($gsvar) || filemtime($gsvar)<$before)
			{
				$skip = false;
			}
		}
		
		//check CNVs
		if (contains($steps, "cn"))
		{
			$cnvs = substr($bam, 0, -4)."_cnvs.seg";
			if (!file_exists($cnvs) || filemtime($cnvs)<$before)
			{
				$skip = false;
			}
		}
		
		if($skip)
		{
			print "$i/$count: Skipped '$ps' because it was analyzed after ".date("m.d.Y", $before)."\n";
			continue;
		}
	}

	//perform analysis
	if ($mode=="print")
	{
		print "$i/$count: Analyzing '$folder'...\n";
	}
	else if ($mode=="default")
	{
		print "$i/$count: Analyzing '$folder'...\n";
	
		$args = array();
		$args[] = "-name $ps";
		$args[] = "-folder $folder";
		$args[] = "-steps $steps";
		$args[] = "-annotation_only"
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
	else if ($mode=="sge")
	{
		print "$i/$count: Queuing '$folder'.\n";
		$args = array();
		$args[] = "-steps $steps";
		$args[] = "-annotation_only";
		list($stdout, $stderr, $return) = $parser->execTool("NGS/db_queue_analysis.php", "-type 'single sample' -samples $ps -args '".implode(" ", $args)."'", false);
		if ($return!=0)
		{
			print "  Error occurred:\n";
			$lines = array_merge($stdout, $stderr);
			foreach($lines as $line)
			{
				print "  ".trim($line)."\n";
			}
		}
		
		++$analyzed;
		if (!is_null($max) && $analyzed>=$max)
		{
			print "Maximum number of SGE jobs reached - stopping\n";
			break;
		}
	}	
}

?>