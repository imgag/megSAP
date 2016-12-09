<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

print "##fileformat=VCFv4.1\n";
print "##INFO=<ID=ID,Number=.,Type=String,Description=\"COSMIC accession(s)\">\n";
print "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO\n";

$in = fopen("php://stdin", "r");
while(!feof($in))
{
	$line =  trim(fgets($in));
	if ($line=="" || $line[0]=="#") continue;
	
	$line = explode("\t", $line);
	$line[7] = "ID=".$line[2];
	$line[2] = ".";
	print implode("\t", $line)."\n";
}


?>