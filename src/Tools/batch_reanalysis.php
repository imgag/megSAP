<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("batch_reanalysis", "Creates commands for batch re-analysis of single samples.");
$parser->addStringArray("samples", "Processed sample names.", false);
$parser->addString("steps", "Analysis steps to perform.", false);
$parser->addInt("threads", "Number of threads used.", true);
$parser->addFlag("annotation_only", "Performs only a reannotation of the already created variant calls.");
$parser->addString("before", "Only samples analyzed before the date are reanalyzed (considers 'steps', format 'DD.MM.YYYY').", true);
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

//init
$analyze_script = "php ".realpath(dirname($_SERVER['SCRIPT_FILENAME'])."/../Pipelines/analyze.php");

foreach($samples as $ps)
{
	$ps = trim($ps);
	
	//determine folder name
	$db_conn = DB::getInstance("NGSD");
	$info = get_processed_sample_info($db_conn, $ps, false);
	if (is_null($info))
	{
		fwrite(STDERR, "Skipped '$ps' because it was not found in NGSD!\n");
		continue;
	}
	$folder = $info['ps_folder'];
	
	//check folder exists
	if(!file_exists($folder))
	{
		fwrite(STDERR, "Skipped '$ps' because folder does not exist: $folder\n");
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
		
		$base = dirname($bam)."/".basename2($bam);
		
		//check VCF/GSvar
		if (contains($steps, "vc"))
		{
			if (!$annotation_only)
			{
				$vcf = $base."_var.vcf.gz";
				if (!file_exists($vcf) || filemtime($vcf)<$before)
				{
					$skip = false;
				}
			}
			
			$gsvar = $base.".GSvar";
			if (!file_exists($gsvar) || filemtime($gsvar)<$before)
			{
				$skip = false;
			}
		}
		
		//check CNVs
		if (contains($steps, "cn"))
		{
			$cnvs = $base."_cnvs_clincnv.tsv";
			if (!file_exists($cnvs) || filemtime($cnvs)<$before)
			{
				$skip = false;
			}
		}
		
		//check SVs
		if (contains($steps, "sv"))
		{
			$svs = $base."_manta_var_structural.bedpe";
			if (!file_exists($svs) || filemtime($svs)<$before)
			{
				$skip = false;
			}
		}
		
		if($skip)
		{
			fwrite(STDERR, "Skipped '$ps' because it was analyzed after ".date("m.d.Y", $before)."\n");
			continue;
		}
	}

	//create command
	$args = array();
	$args[] = "-name $ps";
	$args[] = "-folder $folder";
	$args[] = "-steps $steps";
	if (isset($threads))
	{
		$args[] = "-threads $threads";
	}
	if ($annotation_only)
	{
		$args[] = "-annotation_only";
	}
	
	print "{$analyze_script} ".implode(" ", $args)."\n";
}

?>