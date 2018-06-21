<?php
/** 
	@page find_same_sample_kasp
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("find_same_sample_kasp", "Searches for the sample sample using KASP data.");
$parser->addString("sample", "Sample name.", false);
$parser->addInfile("file", "Sample data file (_converted.tsv).", false);
$parser->addString("folder", "KASP data folder name.", false);
$parser->addInt("max_mm", "Maximum number of mismatches.", true, 2);
$parser->addInt("min_snps", "Minimum number of SNPs.", true, 6);
extract($parser->parse($argv));

//parse SNP set from file name
$parts = explode("_", $file);
foreach($parts as $part)
{
	if (starts_with($part, "Set"))
	{
		$set = $part;
	}
}
print "##Sample: {$sample}\n";
print "##Set: {$set}\n";

//parse sample genotypes from file
$genos = array();
$tmp = file($file);
foreach($tmp as $line)
{
	if (starts_with($line, $sample))
	{
		$genos = explode("\t", trim($line));
		$genos = array_slice($genos, 1);
	}
}
print "##Genotypes of sample: ".implode(" ", $genos)."\n";
if (count($genos)==0)
{
	trigger_error("Found not genotypes for sample '{$sample}' file '{$file}'!", E_USER_ERROR);
}

//search in other files
list($files) = exec2("find $folder -name '*_{$set}_converted.tsv'");
print "##Found ".count($files)." TSV files to check\n";
print "#file\tsample\tmatch\tmismatch\tmatch_perc\n";
foreach($files as $file)
{
	$tmp = file($file);
	foreach($tmp as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		
		$genos2 = explode("\t", $line);
		$sample2 = $genos2[0];
		$genos2 = array_slice($genos2, 1);
		
		$ma = 0;
		$mm = 0;
		for($i=0; $i<count($genos); ++$i)
		{
			if ($genos[$i]=="n/a" || $genos2[$i]=="n/a") continue;

			if ($genos[$i]==$genos2[$i])
			{
				++$ma;
			}
			else
			{
				++$mm;
			}
		}
		if ($ma+$mm<$min_snps) continue;
		if ($mm>$max_mm) continue;
		
		print "$file\t$sample2\t$ma\t$mm\t".number_format(100.0*$ma/($ma+$mm), 2)."\n";
	}
}



?>