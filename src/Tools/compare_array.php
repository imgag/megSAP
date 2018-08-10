<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("compare_array", "Compare SNP genotypes from array and NGS.");
$parser->addInfile("a", "SNP genotype export for CytoScan750K array.", false);
$parser->addInfile("n", "NGS variant list in GSvar format.", false);
$parser->addFlag("mm", "Print mismatches");
extract($parser->parse($argv));

//load chip data
$chip = array();
$file = file($a);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	$parts = explode("\t", $line);
	if ($parts[0]=="Probe Set ID") continue;

	$geno = trim($parts[5]);
	if ($geno=="" || contains($geno, "N")) continue;
	
	$tag = "chr".trim($parts[7]).":".trim($parts[8]);
	$chip[$tag] = $geno;
}
print "Loaded ".count($chip)." genotypes from chip (".basename($a).")\n";

//load ngs data (GSvar)
$ngs = array();
$file = file($n);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	$parts = explode("\t", $line);
	$ref = $parts[3];
	$obs = $parts[4];
	if ($ref=="-" || $obs=="-" || strlen($ref)!=1 || strlen($obs)!=1) continue;
	
	$tag = $parts[0].":".$parts[1];
	$geno = strtr($parts[5], array("wt"=>"{$ref}{$ref}", "het"=>"{$ref}{$obs}", "hom"=>"{$obs}{$obs}"));
	$ngs[$tag] = $geno;
}
print "Loaded ".count($ngs)." genotypes from NGS (".basename($n).")\n";

//compare
$mismatches = array();
$c_overlap = 0;
$c_match = 0;
foreach($chip as $tag => $geno)
{
	if (isset($ngs[$tag]))
	{
		++$c_overlap;
		if ($geno==$ngs[$tag] || $geno[1].$geno[0]==$ngs[$tag])
		{
			++$c_match;
		}
		else
		{
			$mismatches[] = "$tag chip=$geno ngs=".$ngs[$tag];
		}
	}
}
print "Overlapping genotypes: {$c_overlap}\n";
print "Matching genotypes: {$c_match} (".number_format(100.0*$c_match/$c_overlap,2)."%)\n";

//print mismatches
if ($mm)
{
	foreach($mismatches as $mm)
	{
		print "Mismatch: {$mm}\n";
	}
}

?>