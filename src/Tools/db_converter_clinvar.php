<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

function substr_clear($str, $start)
{
	return trim(strtolower(strtr(substr($str, $start), "|", ",")));
}


//locale for iconv below
setlocale(LC_ALL, 'de_DE.UTF8');

//print header
print "##fileformat=VCFv4.1\n";
print "##INFO=<ID=SIG,Number=.,Type=String,Description=\"ClinVar clinical significance\">\n";
print "##INFO=<ID=ACC,Number=.,Type=String,Description=\"ClinVar variation ID\">\n";
print "##INFO=<ID=DISEASE,Number=.,Type=String,Description=\"ClinVar disease annotation\">\n";
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

$debug = false;

//reduce/reformat info
$in = fopen("php://stdin", "r");
while(!feof($in))
{
	$line =  trim(fgets($in));
	if ($line=="" || $line[0]=="#") continue;
	
	$parts = explode("\t", $line);
	
	//skip special chromosomes
	if (chr_check($parts[0], 22, false)===FALSE)
	{
		if ($debug) print "SKIPPED CHR: $line\n";
		continue;
	}
	
	//skip MUT=REF variants
	if ($parts[4]==".")
	{
		if ($debug) print "SKIPPED MUT=REF: $line\n";
		continue;
	}
	
	//extract accession
	$acc = $parts[2];
	
	//extract clinical significance and disease name (both for single variant and variant combinations, e.g. comp-het)
	$infos = explode(";", $parts[7]);
	$sig = "";
	$dis = "";
	$sig_inc = "";
	$dis_inc = "";
	$acc_inc = "";
	foreach($infos as $info)
	{
		if (starts_with($info, "CLNSIG="))
		{
			$sig = substr_clear($info, 7);
		}
		if (starts_with($info, "CLNDN="))
		{
			$dis = substr_clear($info, 6);
		}
		if (starts_with($info, "CLNSIGINCL="))
		{
			list($acc_inc, $sig_inc) = explode(":", substr_clear($info, 11), 2);
		}
		if (starts_with($info, "CLNDNINCL="))
		{
			$dis_inc = substr_clear($info, 10);
		}
	}
	
	$sigs = array();
	$accs = array();
	$diss = array();
	if ($sig!="")
	{
		$sigs[] = $sig;
		$accs[] = $acc;
		$diss[] = $dis;
	}
	if ($sig_inc!="")
	{
		$sigs[] = $sig_inc;
		$accs[] = $acc_inc;
		$diss[] = $dis_inc;
	}
		
	//output
	if  (count($sigs))
	{
		$parts[0] = "chr".$parts[0];
		$parts[2] = ".";
		$parts[5] = ".";
		$parts[6] = ".";
		
		$disease = iconv("utf-8","ascii//TRANSLIT", implode("|", $diss));
		$parts[7] = "SIG=".implode("|", $sigs).";ACC=".implode("|", $accs).";DISEASE=".$disease."";		
		print implode("\t", $parts)."\n";
	}
	else		
	{
		if ($debug) print "SKIPPED NOSIG: $line\n";
	}
}

?>