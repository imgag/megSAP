<?php

include("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

$filter_tags = explode(",", trim($argv[2]));
$min_hits = trim($argv[3]);
$out = $argv[4];

$regs = [];
$h = gzopen2($argv[1], "r");
while(!gzeof($h))
{
	$line = trim(gzgets($h));
	if ($line=="" || $line[0]=="#") continue;
	
	list($chr, $pos, $id, $ref, $obs, $qual, $filter, $info) = explode("\t", $line);
	
	$tag_hit = false;
	foreach($filter_tags as $tag)
	{
		if (contains($filter, $tag)) $tag_hit = true;
	}
	
	if ($tag_hit)
	{
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
$output = array();
foreach($regs as $chr => $tmp)
{
	foreach($tmp as $pos => $count)
	{
		if ($count>=$min_hits)
		{
			$output[] = "$chr\t".($pos-1)."\t".($pos)."\n";
		}
	}
}
file_put_contents($out, $output);

//merge regions
exec2(get_path("ngs-bits")."BedMerge -in $out -out $out");

?>