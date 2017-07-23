<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");


//Fixes several issues with the ClinVar VCF format:
//-splits mutli-allelic variants
//-removed MUT=REF variants

//print header
print "##fileformat=VCFv4.1\n";
print "##INFO=<ID=SIG,Number=.,Type=String,Description=\"ClinVar clinical significance\">\n";
print "##INFO=<ID=ACC,Number=.,Type=String,Description=\"ClinVar accession\">\n";
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

$sig_map = array(255=>"other", 0=>"uncertain_significance", 1=>"not_provided", 2=>"benign", 3=>"likely_benign", 4=>"likely_pathogenic", 5=>"pathogenic", 6=>"drug_response", 7=>"histocompatibility");

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
	
	$infos = explode(";", $line[7]);
	$sigs = null;
	$accs = null;
	$allele_indices = null;
	foreach($infos as $info)
	{
		if (starts_with($info, "CLNSIG="))
		{
			$sigs = explode(",", strtr(substr($info, 7), $sig_map));
		}
		if (starts_with($info, "CLNACC="))
		{
			$accs = explode(",", substr($info, 7));
		}
		if (starts_with($info, "CLNALLE="))
		{
			$allele_indices = explode(",", substr($info, 8));
		}
	}
	
	$alleles = explode(",", $line[3].",".$line[4]);
	
	for($i=0; $i<count($allele_indices); ++$i)
	{
		$allele_i = $allele_indices[$i];
		if ($allele_i==-1 || $allele_i==0) continue; //skip MUT=REF variants
		
		$line[4] = $alleles[$allele_i];
		$line[7] = "SIG=".$sigs[$i].";ACC=".$accs[$i];
		print implode("\t", $line)."\n";
	}	
}



?>