<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//Fixes several issues with the HGMD VCF format:
//-missing 'chr' for chromosomes
//-spaces in the INFO field (otherwise IGV does not load the file)
//-additionally copies the ID field to the INFO field (allows a more compact annotation with VEP --custom argument)
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
	
	//IGV does not support VCF v4.3 > make it compatible with v4.2
	if (starts_with($line, "##fileformat")) $line = "##fileformat=VCFv4.2";
	if (starts_with($line, "##INFO=<ID=PROT,")) $line = "##INFO=<ID=PROT,Number=1,Type=String,Description=\"Protein annotation\">";
	$line = strtr($line, array("%3D"=>"="));
	
	print $line."\n";
}


?>