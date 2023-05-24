<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//Fixes several issues with the HGMD VCF format:
//-missing 'chr' for chromosomes
//-additionally copies the ID field to the INFO field (allows a more compact annotation with VEP --custom argument)
//-removed MUT=REF variants
//-makes the output IGV-compatible

$in = fopen2("php://stdin", "r");
while(!feof($in))
{
	$line =  trim(fgets($in));
	if ($line=="") continue;
	if ($line[0]!="#")
	{
		$line = explode("\t", $line);
		if (contains($line[7], "MUT=REF")) continue;
		
		$line[7] = strtr($line[7], ["%3A"=>":", "%2C"=>"", " "=>"_"]);
		$line[7] = "ID=".$line[2].";".$line[7];
		
		$line = "chr".implode("\t", $line);
	}
	
	//IGV does not support VCF v4.3 > make it compatible with v4.2
	if (starts_with($line, "##fileformat")) $line = "##fileformat=VCFv4.2";
	if (starts_with($line, "##INFO=<ID=PROT,")) $line = "##INFO=<ID=PROT,Number=1,Type=String,Description=\"Protein annotation\">";
	
	//add INFO description for ID column right before header line
	if (starts_with($line, "#CHROM"))
	{
		print "##INFO=<ID=ID,Number=1,Type=String,Description=\"HGMD id\">\n";
	} 
	
	print $line."\n";
}


?>