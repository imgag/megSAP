<?php

include("../../src/Common/all.php");


$fasta = genome_fasta($argv[1]);

//load FAI file
$file = file($fasta.".fai");
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="") continue;
	
	list($chr, $size) = explode("\t", $line);
	
	if (starts_with($chr, "chrGL")) continue;
	if (starts_with($chr, "chrNC")) continue;
	if ($chr == "chrhs37d5") continue;	
	
	//getting chr sequence
	$n_range_start = null;
	list($stdout) = exec2(get_path("samtools")." faidx --length {$size} {$fasta} {$chr}:1-{$size}");
	$seq = trim($stdout[1])."X";
	
	//parsing N regions
	for($i=0; $i<strlen($seq); ++$i)
	{
		if ($seq[$i]=='N')
		{
			if (is_null($n_range_start))
			{
				$n_range_start = $i;
			}
		}
		else
		{
			if (!is_null($n_range_start))
			{
				print "{$chr}\t{$n_range_start}\t{$i}\n";
				$n_range_start = null;
			}
		}
	}
}

?>

