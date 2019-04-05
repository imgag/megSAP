<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("low_cov_regions", "Determines regions that have low coverage in a large part of samples.");
$parser->addInfileArray("in", "Low-coverage BED files for several samples OR one file containing file names.", false);
$parser->addOutfile("out", "Output BED file with low-coverage regions.", false);
//optional
$parser->addInt("percentile", "Percentile of samples with low coverage needed for output.", true, 40);
$parser->addFlag("tumor", "Only process tumor samples, otherwise only non-tumor samples are process).", true, 40);
extract($parser->parse($argv));

//load input file names
if (count($in)==1)
{
	$in = file($in[0]);
}

//process input files
$samples_processed = 0;
$db = DB::getInstance("NGSD");
foreach($in as $bed)
{
	$bed = trim($bed);
	if ($bed=="") continue;
	
	//check quality
	$base = basename($bed, "_lowcov.bed");
	$info = get_processed_sample_info($db, $base, false);
	if(is_null($info))
	{
		print "skipping: $base (not found in NGDS)\n";
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
	
	//filter for female samples (otherwise results for chrX are wrong)
	if ($info['gender']!="female")
	{
		print "skipping: $base (not 'female' gender)\n";
		continue;
	}
	
	++$samples_processed;
	
	//process low-coverage file
	print "processing: $base\n";
	$file = file($bed);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		list($chr, $start, $end) = explode("\t", $line);			
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

//write output BED file
$output = array();
foreach($low as $pos => $count)
{
	if ($count>=$thres)
	{
		list($chr, $pos) = explode("_", $pos);
		$output[] = "$chr\t$pos\t".($pos+1)."\t$count\n";
	}
}
file_put_contents($out, $output);

//merge regions
$parser->exec(get_path("ngs-bits")."BedMerge", "-in $out -out $out", false);

//annotate with gene names
$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $out -clear -out $out", false);

//output
print "input files: ".count($in)."\n";
print "threshold: {$thres} of {$samples_processed} processed samples\n";
list($stdout, $stderr) = exec2(get_path("ngs-bits")."BedInfo -in $out | grep Bases");
list(,$bases) = explode(":", $stdout[0]);
print "low-coverage bases: ".number_format($bases, 0)."\n";
 

?>
