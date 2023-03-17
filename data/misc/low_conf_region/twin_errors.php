<?php

include("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

$out = $argv[1];
$vcfs = array_slice($argv, 2);

$regs = [];
foreach($vcfs as $vcf)
{
	print "processing $vcf...\n";
	$h = gzopen2($vcf, "r");
	while(!gzeof($h))
	{
		$line = trim(gzgets($h));
		if ($line=="" || $line[0]=="#") continue;
		
		list($chr, $pos, $id, $ref, $obs, $qual, $filter, $info) = explode("\t", $line);
		
		$ref_l = strlen($ref);
		$obs_l = strlen($obs);
		$start = $pos;
		if ($ref_l==1 && $obs_l==1) //SNV
		{
			$end = $pos;
		}
		else if ($ref_l==1) //insertions
		{
			$end = $pos+1;
		}
		else //deletion or complex indel
		{
			$end = $pos+$ref_l-1;
		}
		
		for($p=$start; $p<=$end; ++$p)
		{
			if (!isset($regs[$chr][$p]))
			{
				$regs[$chr][$p] = 0;
			}
			$regs[$chr][$p] += 1;
		}
	}
}

//write output BED file
$c_bases = 0;
$output = array();
foreach($regs as $chr => $tmp)
{
	foreach($tmp as $pos => $count)
	{
		if ($count>=2)
		{
			$output[] = "$chr\t".($pos-1)."\t".($pos)."\n";
			++$c_bases;
		}
	}
}
file_put_contents($out, $output);
print "bases: $c_bases\n";

//merge regions
exec2(get_path("ngs-bits")."BedMerge -in $out -out $out");

?>