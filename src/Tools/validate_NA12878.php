<?php
/** 
	@page validate_NA12878 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("validate_NA12878", "\$Rev: 712 $", "Validates that all NA12878 variants are found and that no additional variants are found.");
$parser->addInfile("vcf", "Input variant list of sequencing experiment (VCF.GZ format).", false);
$parser->addInfile("bam", "Mapped reads file (BAM format).", false);
$parser->addInfile("roi", "Input target region file (BED format).", false);
$parser->addFlag("lc", "Also evaluate variants in low-confidence regions");
extract($parser->parse($argv));

//returns the base count of a BED file
function get_bases($filename)
{
	global $parser;
	
	list($stdout) = $parser->exec(get_path("ngs-bits")."BedInfo", "-in $filename", true);
	$hits = array_containing($stdout, "Bases ");
	$parts = explode(":", $hits[0]);
	return trim($parts[1]);
}

//returns the variants of a VCF file in the given ROI
function get_variants($vcf_gz, $roi)
{
	global $parser;
	
	$tmp = $parser->tempFile(".vcf");
	$parser->exec("tabix", "-B $vcf_gz $roi > $tmp", true);
	
	$output = array();
	$file = file($tmp);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;

		$p = explode("\t", $line);
		$var = array();
		$var["QUAL"] = $p[5];
		
		//INFO
		$info = explode(";",$p[7]);
		foreach($info as $entry)
		{
			if (starts_with($entry, "MQM="))
			{
				$var['MQM'] = substr($entry, 4);
			}
		}
		
		//FORMAT, SAMPLE
		$p8 = explode(":", $p[8]);
		$p9 = explode(":", $p[9]);
		$idx = array_search("GT", $p8);
		if ($idx!==FALSE)
		{
			$gt = strtr($p9[$idx], array("|"=>"/", "."=>"0"));
			if ($gt=="1/0") $gt="0/1";
			$var['GT'] = $gt;
		}
		$idx = array_search("DP", $p8);
		if ($idx!==FALSE) $var['DP'] = $p9[$idx];
		$idx = array_search("AO", $p8);
		if ($idx!==FALSE) $var['AO'] = $p9[$idx];
		$idx = array_search("RO", $p8);
		if ($idx!==FALSE) $var['RO'] = $p9[$idx];

		$pos = $p[0].":".$p[1]." ".$p[3].">".$p[4];
		$output[$pos] = $var;
	}
	return $output;
}

//returns the depth at a specified position in the BAM file
function get_depth($pos, $bam)
{
	global $parser;

	list($chr, $start) = explode(" ", strtr($pos, ":", " "));
	
	list($stdout) = $parser->exec(get_path("samtools"), "depth -Q 1 -q 20 -r $chr:$start-$start $bam", true);
	
	if(count($stdout)==0) return 0;
	
	list(, , $depth) = explode("\t", trim($stdout[0]));
	return $depth;
}

//Target region base statistics
print "#Target region     : $roi\n";
$bases = get_bases($roi);
print "#Bases             : $bases\n";
$roi_hc = $parser->tempFile(".bed");
$parser->exec(get_path("ngs-bits")."BedIntersect", "-in $roi -in2 ".get_path("data_folder")."/NA12878/high_conf_regions.bed -out $roi_hc", true);
$based_hc = get_bases($roi_hc);
print "#High-conf bases   : $based_hc (".number_format(100*$based_hc/$bases, 2)."%)\n";
print "#\n";

//get NA12878 variants in ROI
print "#Variant list: $vcf\n";
print "#Notice: Evaluating ".($lc ? "low- and high confidence regions" : "high confidence regions only")."!\n";
$roi_used = $lc ? $roi : $roi_hc;
$found = get_variants($vcf, $roi_used);
print "#Variants observed : ".count($found)."\n";
$expected = get_variants(get_path("data_folder")."/NA12878/high_conf_variants.vcf.gz", $roi_used);
print "#Variants expected : ".count($expected)."\n";

//find missing NA12878 variants or variants with genotype mismatch
$var_diff = array();
foreach($expected as $pos => $var)
{
	if (!isset($found[$pos]))
	{
		$exp = array();
		$exp["DP"] = get_depth($pos, $bam);
		
		$var_diff[$pos] = array("-", $var, $exp);
	}
	else if ($found[$pos]['GT']!=$var['GT'])
	{		
		$var_diff[$pos] = array("g", $var, $found[$pos]);
	}
}

//find unexpected variants
foreach($found as $pos => $var)
{
	if (!isset($expected[$pos]))
	{
		$var_diff[$pos] = array("+", array(), $var);
	}
}

//sort variant list
function pos_gt($a, $b)
{
	list($c_a, $p_a) = explode(" ", strtr($a, array("chr"=>"", ":"=>" ")));
	list($c_b, $p_b) = explode(" ", strtr($b, array("chr"=>"", ":"=>" ")));
	
	if ($c_a==$c_b)
	{
		if ($p_a==$p_b) return 0;
		return ($p_a<$p_b) ? -1 : 1;
	}
	
	return ($c_a<$c_b) ? -1 : 1;
}
uksort($var_diff, "pos_gt");

//print output
$sample = substr(basename($vcf), 0, 10);
function get_prop($props, $name, $digits = -1)
{
	if (!isset($props[$name])) return "";
	if ($digits>=0) number_format($props[$name], $digits);
	return $props[$name];
}
print "#sample\tstatus\tpos\tvariant\texplaination\tref_GT\tobs_GT\tobs_DP\tobs_QUAL\tobs_MQM\tobs_AO\n";
foreach($var_diff as $pos => $data)
{
	list($type, $ref, $obs) = $data;
	$expl = array();
	if ($type=="-")
	{
		if(get_prop($obs, 'DP')<3)
		{
			$expl[] = "DP<3";
		}
		else if(get_prop($obs, 'DP')<10)
		{
			$expl[] = "DP<10";
		}
	}
	if ($type=="g")
	{
		if(get_prop($obs, 'DP')<10)
		{
			$expl[] = "DP<10";
		}
	}
	if ($type=="+")
	{
		if (get_prop($obs, 'MQM')<50)
		{	
			$expl[] = "MQM<50";
		}
	}
	print "$sample\t$type\t".strtr($pos, " ", "\t")."\t".implode(", ", $expl)."\t".get_prop($ref, "GT")."\t".get_prop($obs, "GT")."\t".get_prop($obs, "DP")."\t".get_prop($obs, "QUAL", 2)."\t".get_prop($obs, "MQM", 2)."\t".get_prop($obs, "AO")."\n";
}


?>