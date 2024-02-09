<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

function substr_clear($str, $start, $use_tab_as_sep=false)
{
	return trim(strtr(substr($str, $start), array("|"=> ($use_tab_as_sep ? "\t" : ","))));
}

function cons2value($cons, $line)
{
	//in ClinVar 20220328 there are extra underscores...
	$cons = trim($cons, "_");
	
	//remove consequence modifiers
	$cons = strtr($cons, [",_low_penetrance"=>""]);

	if ($cons=="pathogenic") return 5;
	if ($cons=="likely_pathogenic") return 4;
	if ($cons=="uncertain_significance") return 3;
	if ($cons=="likely_benign") return 2;
	if ($cons=="benign") return 1;
	if ($cons=="likely_risk_allele") return -1;
	if ($cons=="uncertain_risk_allele") return -1;
	if ($cons=="established_risk_allele") return -1;
	trigger_error("Unknown consequence '$cons' in function 'cons2value' in line:\n$line", E_USER_ERROR);
}

function value2cons($value)
{
	if ($value<1 || $value>5)
	{
		trigger_error("Invalid consequence value '$value' in function 'value2cons'!", E_USER_ERROR);
	}
	if ($value>4.5) return "Pathogenic";
	if ($value>3.5) return "Likely_pathogenic";
	if ($value<1.5) return "Benign";
	if ($value<2.5) return "Likely_benign";
	return "Uncertain_significance";
}
//locale for iconv below
setlocale(LC_ALL, 'de_DE.UTF8');

//print header
print "##fileformat=VCFv4.1\n";
print "##INFO=<ID=DETAILS,Number=.,Type=String,Description=\"ClinVar disease/significance annotation\">\n";
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

$debug = false;

//reduce/reformat info
$in = fopen2("php://stdin", "r");
while(!feof($in))
{
	$line =  trim(fgets($in));
	if ($line=="" || $line[0]=="#") continue;
	
	$parts = explode("\t", $line);
	list($chr, $pos, $acc, $ref, $alt) = $parts;
	$chr = "chr".$chr;
	
	//skip special chromosomes
	if (chr_check($chr, 22, false)===FALSE)
	{
		if ($debug) print "SKIPPED CHR: $line\n";
		continue;
	}
	
	//skip ALT=REF variants
	if ($alt==".")
	{
		if ($debug) print "SKIPPED MUT=REF: $line\n";
		continue;
	}
	
	//skip invalid alternative alleles
	if (!preg_match("/^[NACGT]+$/", $alt))
	{
		if ($debug) print "SKIPPED invalid alternative allele: $line\n";
		continue;
	}
	
	//prepend padding base for complex variants
	if ((strlen($ref)>1 || strlen($alt)>1) && strlen($ref)!=strlen($alt) && $ref[0]!=$alt[0])
	{
		$pos -= 1;
		$base = get_ref_seq("GRCh38", $chr, $pos, $pos);
		$ref = $base.$ref;
		$alt = $base.$alt;
	}
	
	//extract clinical significance and disease name (both for single variant and variant combinations, e.g. comp-het)
	$infos = explode(";", $parts[7]);
	$sig = "";
	$sig_conf = "";
	$replace_msg = "";
	$dis = "";
	$sig_acc_inc = [];
	$dis_inc = "";
	$genes = [];
	foreach($infos as $info)
	{
		if (starts_with($info, "CLNSIG="))
		{
			$sig = substr_clear($info, 7);
		}
		if (starts_with($info, "CLNSIGCONF="))
		{
			$sig_conf = strtolower(substr_clear($info, 11, true));
		}
		if (starts_with($info, "CLNDN="))
		{
			$dis = substr_clear($info, 6);
		}
		if (starts_with($info, "CLNSIGINCL="))
		{
			$parts2 = explode(",", substr_clear($info, 11));
			foreach($parts2 as $part)
			{
				$sig_acc_inc[] = explode(":", $part.":");
			}
		}
		if (starts_with($info, "CLNDNINCL="))
		{
			$dis_inc = substr_clear($info, 10);
		}
		if (starts_with($info, "GENEINFO="))
		{
			$tmp = explode("|",substr_clear($info, 9));
			foreach($tmp as $entry)
			{
				$gene = trim(explode(":", $entry)[0]);
				if (!in_array($gene, $genes)) $genes[] = $gene;
			}
		}
	}
	
	//handle conflicting interpretations that are not conflicting
	if ($sig_conf!="")
	{
		if (!(contains($sig_conf, "benign") && contains($sig_conf, "pathogenic")))
		{
			//calcualate the mean consequence
			$cons_counts = array(1=>0, 2=>0, 3=>0, 4=>0, 5=>0);
			$sig_conf = explode("\t", $sig_conf);
			foreach($sig_conf as $entry)
			{
				$entry = trim($entry);
				$entry = substr($entry, 0, -1);
				list($cons, $count) = explode("(", $entry);
				$value = cons2value($cons, $line);
				if ($value==-1) continue; //skip risk variants, etc.
				$cons_counts[$value] += $count;
			}
			
			$counts = [];
			$avg = 0;
			foreach($cons_counts as $cons=>$count)
			{
				$avg += $cons * $count;
				if ($count!=0)
				{
					$counts[] = "{$count}x{$cons}";
				}
			}
			$avg /= array_sum(array_values($cons_counts));
			$sig = value2cons($avg);
			$replace_msg = "(replaced_by_megSAP_based_on_classifications:".implode(",", $counts).")";
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
	foreach($sig_acc_inc as list($acc_inc, $sig_inc))
	{
		$sigs[] = $sig_inc;
		$accs[] = $acc_inc;
		$diss[] = $dis_inc;
	}
	
	//skip variant without significance
	if (count($sigs)==0)
	{
		if ($debug) print "SKIPPED NOSIG: $line\n";
		continue;
	}

	//output
	for($i=0; $i<count($sigs); ++$i)
	{
		$info = "DETAILS=".vcf_encode_url_string(strtolower($sigs[$i]).($i==0 ? $replace_msg : "")." DISEASE=".iconv("utf-8", "ascii//TRANSLIT", $diss[$i])." GENE=".implode(",", $genes));
		print implode("\t", [$chr, $pos, $accs[$i], $ref, $alt, ".", ".", $info])."\n";
	}
}

?>