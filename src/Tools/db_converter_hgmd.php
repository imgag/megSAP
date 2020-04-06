<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//Fixes several issues with the HGMD VCF format:
//-missing 'chr' for chromosomes
//-spaces in the INFO field (otherwise IGV does not load the file)
//-additionally copies the ID field to the INFO field (allows a more compact annotation with VEP --custom argument)
//-removed MUT=REF variants

$in = fopen2("php://stdin", "r");
while(!feof($in))
{
	$line =  trim(fgets($in));
	if ($line=="") continue;
	if ($line[0]!="#")
	{
		$line = explode("\t", $line);
		if (contains($line[7], "MUT=REF")) continue;
		//replace invalid characters in INFO column with URL encoded string
		$info_entries = array();
		$info_entries[] = "ID=".$line[2];
		foreach (explode(";", $line[7]) as $key_value_pair) 
		{
			list($key, $value) = explode("=", $key_value_pair, 2);
			$info_entries[] = $key."=".vcf_encode_url_string($value);
		}
		$line[7] = implode(";", $info_entries);
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