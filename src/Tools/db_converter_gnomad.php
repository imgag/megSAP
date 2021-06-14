<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$handle = fopen2("php://stdin", "r");
while($line = fgets($handle))
{
	$line = trim($line);
	if ($line=="") continue;
	if($line[0]=="#")
	{
		if (in_array("-header", $argv))
		{
			if($line[1]=="#")
			{
				//comment section

				//skip INFO ID headers (not present in output file)
				if (starts_with($line, "##INFO=<ID=")) continue;

				//print comment line
				print $line."\n";
			}
			else
			{
				//actual header line

				// add INFO column header:
				print "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in samples\">\n";
				print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency in samples\">\n";
				print "##INFO=<ID=Hom,Number=A,Type=Integer,Description=\"Count of homozygous individuals in samples\">\n";
				print "##INFO=<ID=Hemi,Number=A,Type=Integer,Description=\"Alternate allele count for male samples\">\n";

				//print header line
				print $line."\n";
			}
			
		}
		continue;
	}
	
	$parts = explode("\t", $line);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = $parts;
	$info = explode(";", $info);
	
	//get required infos
	$is_chrx = $chr=="chrX";
	$nonpar = in_array("nonpar", $info);
	$af = null;
	$an = null;
	$hom = null;
	$hemi = null;
	foreach ($info as $entry)
	{
		if (starts_with($entry, "AF="))
		{
			$af = substr($entry, 3);
		}
		else if (starts_with($entry, "AN="))
		{
			$an = substr($entry, 3);
		}
		else if (starts_with($entry, "nhomalt="))
		{
			$hom = substr($entry, 8);
		}
		else if ($is_chrx && $nonpar && starts_with($entry, "AC_male="))
		{
			$hemi = substr($entry, 8);
		}
	}
	
	$info_new = array();
	$info_new[] = "AN=".$an;
	$info_new[] = "Hemi=".($is_chrx && $nonpar ? $hemi : ".");
	$info_new[] = "Hom=".$hom;
	$info_new[] = "AF=".($an<200 ? '0.0' : number_format($af, 5));

	print "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t".implode(";", $info_new)."\n";
}
fclose($handle);

?>