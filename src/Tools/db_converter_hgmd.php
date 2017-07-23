<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//Fixes several issues with the HGMD VCF format:
//-missing 'chr' for chromosomes
//-spaces in the INFO field (otherwise IGV/SnpSift does not load the file)
//-additionally copies the ID field to the INFO field (needed for SnpSift annotation tool)
//-removed MUT=REF variants

$in = fopen("php://stdin", "r");
while(!feof($in))
{
	$line =  trim(fgets($in));
	if ($line=="") continue;
	if ($line[0]!="#")
	{
		$line = explode("\t", $line);
		if (contains($line[7], "MUT=REF")) continue;
		$line[7] = "ID=".$line[2].";".strtr($line[7], array(" "=>"_"));
		$line = "chr".implode("\t", $line);
	}
	print $line."\n";
}


?>