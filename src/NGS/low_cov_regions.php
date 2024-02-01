<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("low_cov_regions", "Determines regions that have low coverage in a large part of samples.");
$parser->addInfileArray("in", "Low-coverage BED files for several samples OR one file containing file names.", false);
$parser->addOutfile("out", "Output BED file with low-coverage regions.", false);
//optional
$parser->addInt("percentile", "Percentile of samples with low coverage needed for output.", true, 40);
$parser->addFlag("tumor", "If set, only tumor samples are processed. Otherwise only non-tumor samples are processed.");
extract($parser->parse($argv));

//load input file names
if (count($in)==1)
{
	$in = file($in[0]);
}

//process input files
$samples_processed = 0;
$samples_processed_x = 0;
$samples_processed_y = 0;
$db = DB::getInstance("NGSD");
foreach($in as $bed)
{
	$bed = trim($bed);
	if ($bed=="") continue;
	
	//check quality
	$base = basename($bed, ".bed");
	$info = get_processed_sample_info($db, $base, false);
	if(is_null($info))
	{
		print "skipping: $base (not found in NGSD)\n";
		continue;
	}
	if ($info['ps_quality']!="good")
	{
		print "skipping: $base (not 'good' quality)\n";
		continue;
	}
	
	//filter out or select tumor samples
	if ($tumor && !$info['is_tumor'])
	{
		print "skipping: $base (is non-tumor sample)\n";
		continue;
	}
	else if (!$tumor && $info['is_tumor'])
	{
		print "skipping: $base (is tumor sample)\n";
		continue;
	}
	
	//count number of samples
	$gender = $info['gender'];
	++$samples_processed;
	if ($gender=="female") ++$samples_processed_x;
	if ($gender=="male") ++$samples_processed_y;
	
	//process low-coverage file
	print "processing: $base\n";
	$file = file($bed);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		list($chr, $start, $end) = explode("\t", $line);
		
		if ($chr=="chrX" && $gender!="female") continue;
		if ($chr=="chrY" && $gender!="male") continue;
		
		for($p=$start; $p<$end; ++$p)
		{
			if (!isset($low[$chr."_".$p]))
			{
				$low[$chr."_".$p] = 0;
			}
			$low[$chr."_".$p] += 1;
		}
	}
}

//calculate threshold
$thres = number_format($samples_processed * $percentile / 100.0, 2);
$thres_x = number_format($samples_processed_x * $percentile / 100.0, 2);
$thres_y = number_format($samples_processed_y * $percentile / 100.0, 2);

//write output BED file
$output = [];
foreach($low as $pos => $count)
{
	list($chr, $pos) = explode("_", $pos);
	if($chr=="chrX" && $count>=$thres_x)
	{
		$output[] = "$chr\t$pos\t".($pos+1)."\n";
	}
	else if($chr=="chrY" && $count>=$thres_y)
	{
		$output[] = "$chr\t$pos\t".($pos+1)."\n";
	}
	else if ($count>=$thres)
	{
		$output[] = "$chr\t$pos\t".($pos+1)."\n";
	}
}
file_put_contents($out, $output);

//pre-merge BED file (otherwise it can be too large for BedMerge)
$i_chr = "";
$i_start = -1;
$i_end = -1;
$output = [];
foreach(file($out) as $line)
{
	$line = trim($line);
	if ($line=="") continue;
	
	list($chr, $start, $end) = explode("\t", $line);
	if ($i_start==-1) //init interval
	{
		$i_chr = $chr;
		$i_start = $start;
		$i_end = $end;
	}
	else
	{
		if ($i_chr==$chr && $i_end==$start) //next base > extend interval
		{
			$i_end = $end;
		}
		else //end of interval > write output and init next interval
		{
			$output[] = "$i_chr\t$i_start\t$i_end\n";
			
			$i_chr = $chr;
			$i_start = $start;
			$i_end = $end;
		}
	}
}
$output[] = "$i_chr\t$i_start\t$i_end\n";
file_put_contents($out, $output);


//merge regions
$parser->exec(get_path("ngs-bits")."BedMerge", "-in $out -out $out", false);

//annotate with gene names
$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $out -clear -out $out", false);

//output
print "input files: ".count($in)."\n";
print "threshold: {$thres} of {$samples_processed} processed samples\n";
print "threshold chrX: {$thres_x} of {$samples_processed_x} processed females\n";
print "threshold chrY: {$thres_y} of {$samples_processed_y} processed males\n";
list($stdout, $stderr) = exec2(get_path("ngs-bits")."BedInfo -in $out | grep Bases");
list(,$bases) = explode(":", $stdout[0]);
print "low-coverage bases: ".number_format($bases, 0)."\n";
 

?>
