<?php

include("../../src/Common/all.php");

//load data
function load($filename, &$chrs, &$breaks)
{
	$output = array();
	$h = fopen2($filename, "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || starts_with($line, "Chromosome")) continue;
		
		list($chr, $start, $end, $hist, $af) = explode("\t", $line);
		
		$output[$chr][] = array($start, $end, $hist, $af);
		
		if (!in_array($chr, $chrs)) $chrs[] = $chr;
		$breaks[$chr][] = $start;
		$breaks[$chr][] = $end;
	}
	return $output;
}
$chrs = array();
$breaks = array();
$f1 = load($argv[1], $chrs, $breaks);
$f2 = load($argv[2], $chrs, $breaks);
$source2 = basename($argv[2], ".igv");

print "Chromosome	Start	End	CN histogram (0-10)	AF TruSeqPCRfree\n";

//process
foreach($chrs as $chr)
{
	//determine borders
	$borders = $breaks[$chr];
	sort($borders);
	$borders = array_unique($borders);
	$borders = array_values($borders);
	
	//for each range, determine and print max
	for ($i=0; $i<count($borders)-1; ++$i)
	{
		$start = $borders[$i];
		$end = $borders[$i+1];
		$max_af = 0.0;
		$max_match = null;
		
		if (isset($f1[$chr]))
		{
			foreach($f1[$chr] as list($start_f, $end_f, $hist, $af))
			{
				$inter = range_intersect($start, $end, $start_f, $end_f);
				if ($inter!==FALSE && $inter[1]-$inter[0]>0)
				{
					if ($af>$max_af)
					{
						$max_match = array($hist, $af, false);
					}
					//print "$start, $end - match ($start_f-$end_f) - $af\n";
				}
			}
		}
		if (isset($f2[$chr]))
		{
			foreach($f2[$chr] as list($start_f, $end_f, $hist, $af))
			{
				$inter = range_intersect($start, $end, $start_f, $end_f);
				if ($inter!==FALSE && $inter[1]-$inter[0]>0)
				{
					if ($af>$max_af)
					{
						$max_match = array($hist, $af, true);
					}
					//print "$start, $end - match ($start_f-$end_f) - $source2 $af\n";
				}
			}
		}
		
		if (!is_null($max_match))
		{
			$hist = strtr(trim($max_match[0]), " ", ",");
			$af = trim($max_match[1]);
			print "{$chr}\t{$start}\t{$end}\t{$hist}".($max_match[2] ? " (CNPs)" : "")."\t{$af}\n";
		}
	}
}

?>