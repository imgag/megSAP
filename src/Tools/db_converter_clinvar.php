<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//print header
print "##fileformat=VCFv4.1\n";
print "##INFO=<ID=SIG,Number=.,Type=String,Description=\"ClinVar clinical significance\">\n";
print "##INFO=<ID=ACC,Number=.,Type=String,Description=\"ClinVar accession\">\n";
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

$sig_map = array(255=>"other", 0=>"uncertain_significance", 1=>"not_provided", 2=>"benign", 3=>"likely_benign", 4=>"likely_pathogenic", 5=>"pathogenic", 6=>"drug_response", 7=>"histocompatibility", ","=>"|");

//reduce/reformat  info
$in = fopen("php://stdin", "r");
while(!feof($in))
{
	$line =  trim(fgets($in));
	if ($line=="" || $line[0]=="#") continue;
	
	$line = explode("\t", $line);
	$line[0] = "chr".$line[0];
	$line[2] = ".";
	$line[5] = ".";
	$line[6] = ".";
	
	$annos = array();
	$infos = explode(";", $line[7]);
	foreach($infos as $info)
	{
		if (starts_with($info, "CLNSIG="))
		{
			$annos[] = "SIG=".strtr(substr($info, 7), $sig_map);
		}
		if (starts_with($info, "CLNACC="))
		{
			$annos[] = "ACC=".strtr(substr($info, 7), array(","=>"|"));
		}
	}
	$line[7] = implode(";", $annos);
	
	print implode("\t", $line)."\n";
}


?>