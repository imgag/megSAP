<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("low_cov_regions", "Determines regions that have low coverage in a large part of samples.");
$parser->addInfileArray("in", "Low-coverage BED files for several samples OR one file containing file names.", false);
$parser->addOutfile("out", "Output BED file with low-coverage regions.", false);
//optional
$parser->addInt("percentile", "Percentile of samples with low coverage needed for output.", true, 40);
$parser->addFlag("tumor", "Only process tumor samples, otherwise only non-tumor samples are process).");
$parser->addFlag("only_females", "Only process female samples (otherwise gaps for chrX are over-estimated).");
extract($parser->parse($argv));

//load input file names
if (count($in)==1)
{
	$in = file($in[0]);
}

//process input files
$samples_to_process = 0;
$db = DB::getInstance("NGSD");

$bed_files_filtered = array();

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
	if ($only_females && $info['gender']!="female")
	{
		print "skipping: $base (not 'female' gender)\n";
		continue;
	}
	
	// add to BED file list
	$bed_files_filtered[] = $bed;
	++$samples_to_process;
}

//calculate threshold
$thres = number_format($samples_to_process * $percentile / 100.0, 2);

// get chr list
$chromosomes = array('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrY','chrX','chrMT');

// create empty file
file_put_contents($out, array());
foreach ($chromosomes as $ref_chr) 
{
	print "processing: $ref_chr\n";
	$low = array();
	foreach ($bed_files_filtered as $bed) 
	{
		//process low-coverage file
		print "\t |--$bed\n";
		$file = file($bed);
		foreach($file as $line)
		{
			$line = trim($line);
			if ($line=="" || $line[0]=="#") continue;
			list($chr, $start, $end) = explode("\t", $line);
			
			//skip every other chr
			if ($chr != $ref_chr) continue;

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
	file_put_contents($out, $output, FILE_APPEND);
	//merge after each chr to keep the file small
	$parser->exec(get_path("ngs-bits")."BedMerge", "-in $out -out $out", false);
}



//merge regions
$parser->exec(get_path("ngs-bits")."BedMerge", "-in $out -out $out", false);

//annotate with gene names
$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $out -out $out", false);

//output
print "input files: ".count($in)."\n";
print "threshold: {$thres} of {$samples_to_process} processed samples\n";
list($stdout, $stderr) = exec2(get_path("ngs-bits")."BedInfo -in $out | grep Bases");
list(,$bases) = explode(":", $stdout[0]);
print "low-coverage bases: ".number_format($bases, 0)."\n";
 

?>
