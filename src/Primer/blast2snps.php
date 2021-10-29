<?php

/** 
	@page blast2snps
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("blast2snps", "Looks up SNPs in primer pairs.");
$parser->addInfile("in",  "Input primer TXT file.", false);
$parser->addInfile("blast",  "Input blast output file.", false);
$parser->addOutfile("out",  "Output TXT file.", false);
$parser->addString("db_list",  "Comma-separated list of variant DBs to use.", false);

//optional
$parser->addInt("max_len", "Maximum product length.", true, 2000);
$parser->addFloat("min_freq", "Minmum allele frequency of SNPs.", true, 0.001);
$parser->addOutfile("bed", "Creates an optional BED file with the primer positions.", true);
extract($parser->parse($argv));

//function to calculate hit length
function length($h1, $h2)
{
	return max($h1[3],$h1[4],$h2[3],$h2[4])-min($h1[3],$h1[4],$h2[3],$h2[4])+1;
}

//filters entry according to minimum frequency. Returns 'false' if frequency too low.
function filter_by_af($line, $min_freq)
{
	$parts = explode("\t", $line);
	if (count($parts)<8) return false;

	//skp CNVs
	if (starts_with($parts[4], "<CN")) return false;	
	
	$af = 0.0;
	$info = explode(";", $parts[7]);
	foreach($info as $entry)
	{
		if (starts_with($entry, "AF="))
		{
			$af = substr($entry, 3); 
		}
	}
	if ($af<$min_freq) return false;
	$af = number_format(max(0.0001, $af), 4);
	
	$id = $parts[2];
	if ($id==".") $id = "";
	
	$chr = $parts[0];
	if (!starts_with($chr, "chr")) $chr = "chr".$chr;
	
	return "{$chr}\t{$parts[1]}\t{$parts[3]}\t{$parts[4]}\t{$af}\t{$id}";
}

/*
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
TNNI3d_1R	chr19	100.00	24	0	0	1	24	55668736	55668759	7e-05	48.1
TNNI3d_1R	chr19	100.00	18	0	0	3	20	2194693	2194710	0.28	36.2
TNNI3d_1R	chr19	100.00	17	0	0	3	19	5251756	5251772	1.1	34.2
TNNI3d_1R	chr19	100.00	16	0	0	1	16	3820802	3820787	4.4	32.2
TNNI3d_1R	chr19	100.00	16	0	0	4	19	17166842	17166857	4.4	32.2
*/
//load blast results
$file = file($blast);
$blast = array();
foreach($file as $line)
{
	if (starts_with($line, "#")) continue;
	$parts = explode("\t", $line);
	if (count($parts)<11) continue;
	
	list($name, $chr, $identity, $length, , , , , $start, $end, $evalue) = $parts;
	if (!isset($blast[$name])) $blast[$name] = array();
	
	$blast[$name][] = array($chr, $identity, $length, $start, $end, $evalue);
}

//determine DBs to use
$dbs = array();
$db_list = explode(",", $db_list);
foreach($db_list as $db)
{
	if ($db=="dbSNP") $dbs[$db] = get_path("GRCh37_data_folder")."/dbs/1000G/1000g_v5b.vcf.gz";
	else if ($db=="ESP6500") $dbs[$db] = get_path("GRCh37_data_folder")."/dbs/ESP6500/ESP6500SI_EA_SSA137.vcf.gz";
	else if ($db=="ExAC") $dbs[$db] = get_path("GRCh37_data_folder")."/dbs/ExAC/ExAC_r0.3.1.vcf.gz";
	else trigger_error("Unknown database '$db'!", E_USER_ERROR);
}
//print DB versions
foreach($dbs as $db_name => $db_file)
{
	print "$db_name version: ".basename($db_file, ".vcf.gz")."\n";
}

//find best pair matches
$input = file($in);
$output = array();
$output_bed = array();
foreach($input as $line)
{
	if (trim($line)=="" || $line[0]=="#" || starts_with($line, "Name FOR")) continue;
	list($n1, $n2, $s1, $s2, $chr_exp, $len_exp) = explode("\t", trim($line));
	$n1 = trim($n1);
	$n2 = trim($n2);
	$s1 = trim($s1);
	$s2 = trim($s2);
	$chr_exp = chr_trim($chr_exp);
	
	//init output
	$output[] = "################################### PAIR ###################################";
	$output[] = "INPUT: ".trim($line);
	
	$hits1 = array();
	if (isset($blast[$n1]))
	{
		$hits1 = $blast[$n1];
	}
	else
	{
		print "WARNING: no hits found for '$n1'. Sequence '$s1' might be too short.\n";
	}
	$hits2 = array();
	if (isset($blast[$n2]))
	{
		$hits2 = $blast[$n2];
	}
	else
	{
		print "WARNING: no hits found for '$n2'. Sequence '$s2' might be  too short.\n";
	}
	
	//reduce to possible hits (correct chr, length<max_len)
	$possible_hits = array();
	foreach ($hits1 as $h1)
	{
		if (chr_trim($h1[0])!=$chr_exp) continue;
		
		foreach ($hits2 as $h2)
		{
			if (chr_trim($h2[0])!=$chr_exp) continue;
			
			if ($h1[0]==$h2[0])
			{
				$dist = length($h1, $h2);
				if ($dist<$max_len)
				{
					$possible_hits[] = array($h1, $h2);
				}
			}
		}
	}

	if (count($possible_hits)==0)
	{
		$output[] = "*ERROR: No blast hit pair found (correct chr and maximum product length of '$max_len'):";
		foreach($hits1 as $h)
		{
			$output[] = "          H1: ".implode(" ", $h);
		}
		foreach($hits2 as $h)
		{
			$output[] = "          H2: ".implode(" ", $h);
		}
		continue;
	}
	else if (count($possible_hits)>=2)
	{
		$output[] = "*WARNING: Several possible blast hit pairs found. Taking the pair with matching product length and lowest e-value sum!";
		$min_sum = 99999999;
		$min_index = -1;
		for ($i=0; $i<count($possible_hits); ++$i)
		{
			$h1 = $possible_hits[$i][0];
			$h2 = $possible_hits[$i][1];
			$length = length($h1, $h2);
			$evalue = $h1[5] + $h2[5];
			if ($length==$len_exp && $evalue<$min_sum)
			{
				$min_index = $i;
				$min_sum = $evalue;
			}
			
			$output[] = "          $i) hit for: ".implode(" ", $possible_hits[$i][0]);
			$output[] = "          $i) hit rev: ".implode(" ", $possible_hits[$i][1]);
			$output[] = "          $i) e-value: $evalue";
			$output[] = "          $i) length : $length";
		}
		
		if ($min_index==-1)
		{
			$output[] = "*ERROR: No matching pair found!";
			continue;
		}
		
		$possible_hits = array($possible_hits[$min_index]);
	}

	//output primer 1
	$hit1 = $possible_hits[0][0];
	$chr1 = chr_trim($hit1[0]);
	$start1 = min($hit1[3],$hit1[4]);
	$end1 = max($hit1[3],$hit1[4]);
	$output[] = "BLAST: $n1 - chr$chr1:$start1-$end1 perc_aligned:".$hit1[1]." length:".$hit1[2]."/".strlen($s1)." e-value:".$hit1[5];
	
	foreach($dbs as $db_name => $db_file)
	{
		$chr_prefix = ($db_name=="ESP6500" ? "chr" : "");
		list($snps) = $parser->exec("tabix", "$db_file {$chr_prefix}$chr1:$start1-$end1", false);
		foreach($snps as $snp)
		{
			$snp = filter_by_af($snp, $min_freq);
			if ($snp!==FALSE) $output[] = "*SNP $db_name: ".htmlspecialchars($snp);
		}
	}
	
	//output primer 2
	$hit2 = $possible_hits[0][1];
	$chr2 = chr_trim($hit2[0]);
	$start2 = min($hit2[3],$hit2[4]);
	$end2 = max($hit2[3],$hit2[4]);
	$output[] = "BLAST: $n2 - chr$chr2:$start2-$end2 perc_aligned:".$hit2[1]." length:".$hit2[2]."/".strlen($s2)." e-value:".$hit2[5];

	foreach($dbs as $db_name => $db_file)
	{
		$chr_prefix = ($db_name=="ESP6500" ? "chr" : "");
		list($snps) = $parser->exec("tabix", "$db_file {$chr_prefix}$chr2:$start2-$end2", false);
		foreach($snps as $snp)
		{
			$snp = filter_by_af($snp, $min_freq);
			if ($snp!==FALSE) $output[] = "*SNP $db_name: ".htmlspecialchars($snp);
		}
	}
	
	//check given chr/dist
	$len = length($hit1, $hit2);
	$len_dev = abs($len_exp-$len);
	if ($len_dev>5)
	{
		$output[] = "*WARNING: Primer product has wrong length. Should be '$len_exp' is '$len', difference '$len_dev'!";
	}
	
	//create optional BED file
	if (isset($bed))
	{
		$poss = array($start1, $start2, $end1, $end2);
		sort($poss, SORT_NUMERIC);
		$output_bed[] = "chr$chr1\t".$poss[0]."\t".$poss[3]."\t$n1-$n2\t700\t*\t".$poss[1]."\t".$poss[2];
	}
}

file_put_contents($out, implode("\n", $output)."\n");

//write optional BED file
if (isset($bed))
{
	file_put_contents($bed, implode("\n", $output_bed)."\n");
}
?>
