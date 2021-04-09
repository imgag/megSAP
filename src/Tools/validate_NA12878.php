<?php
/** 
	@page validate_NA12878
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("validate_NA12878", "Validates the performance of a sequencing experiment on the GiaB reference sample NA12878.");
$parser->addInfile("vcf", "Input variant list of sequencing experiment (VCF.GZ format).", false);
$parser->addInfile("bam", "Mapped reads of sequencing experiment (BAM format).", false);
$parser->addInfile("roi", "Target region of sequencing experiment (BED format).", false);
$parser->addOutfile("stats", "Append statistics to this file.", false);
//optional
$parser->addString("name", "Name used in the 'stats' output. If unset, the 'vcf' file base name is used.", true);
$parser->addInt("min_dp", "If set, only regions in the 'roi' with at least the given depth are evaluated.", true, 0);
$parser->addInt("max_indel", "Maximum indel size (larger indels are ignored).", true, 0);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addString("ref_sample", "Reference sample to use for validation.", true, "NA12878");
$parser->addFlag("matches", "Do not only show variants that were missed (-), are novel (+) or with genotype mismatch (g), but also show matches (=) in output.");
extract($parser->parse($argv));

//returns the base count of a BED file
function get_bases($filename)
{
	global $parser;
	global $ngsbits;
	
	list($stdout) = $parser->exec("{$ngsbits}BedInfo", "-in $filename", true);
	$hits = array_containing($stdout, "Bases ");
	$parts = explode(":", $hits[0]);
	return trim($parts[1]);
}

//returns the variants of a VCF file in the given ROI
function get_variants($vcf_gz, $roi, $normalize, $max_indel)
{
	global $parser;
	global $ngsbits;
	global $vcflib;
	global $genome;
	
	//get variants
	$tmp = $parser->tempFile(".vcf");
	$pipeline = [];
	$pipeline[] = array("zcat", $vcf_gz);
	$pipeline[] = array("{$ngsbits}VcfFilter", "-reg {$roi}");
	if ($normalize)
	{
		$pipeline[] = array("{$vcflib}vcfbreakmulti", "");
		$pipeline[] = array("{$ngsbits}VcfLeftNormalize", "-ref $genome");
	}
	$pipeline[] = array("{$ngsbits}VcfStreamSort", "-out {$tmp}");
	$parser->execPipeline($pipeline, "variant extraction");
	
	//put together output
	$output = array();	
	$file = file($tmp);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;

		//get variant infos
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = explode("\t", $line);
		if (!starts_with($chr, "chr")) $chr = "chr".$chr;
		
		//compile output
		$var = array();
		$var["QUAL"] = $qual;
		$info = explode(";",$info);
		foreach($info as $entry)
		{
			if (starts_with($entry, "MQM="))
			{
				$var['MQM'] = substr($entry, 4);
			}
		}
		
		$format = explode(":", $format);
		$sample = explode(":", $sample);
		$idx = array_search("GT", $format);
		if ($idx!==FALSE)
		{
			$gt = strtr($sample[$idx], array("|"=>"/", "."=>"0"));
			if ($gt=="1/0") $gt="0/1";
			$var['GT'] = $gt;
		}
		$idx = array_search("DP", $format);
		if ($idx!==FALSE) $var['DP'] = $sample[$idx];
		$idx = array_search("AO", $format);
		if ($idx!==FALSE) $var['AO'] = $sample[$idx];
		$idx = array_search("RO", $format);
		if ($idx!==FALSE) $var['RO'] = $sample[$idx];
		$idx = array_search("GQ", $format);
		if ($idx!==FALSE) $var['GQ'] = $sample[$idx];
		if (strlen($ref)>1 || strlen($alt)>1)
		{
			$var["TYPE"] = "INDELS";
			if ($max_indel>0 && (strlen($ref)>$max_indel || strlen($alt)>$max_indel)) continue;
		}
		else
		{
			$var["TYPE"] = "SNVS";
		}
		
		$tag = "{$chr}:{$pos} {$ref}>{$alt}";
		if (isset($output[$tag])) trigger_error("Variant $tag twice in $vcf_gz!", E_USER_ERROR);
		$output[$tag] = $var;
	}
	
	return $output;
}

//returns the depth at a specified position in the BAM file
function get_depth($pos, $bam)
{
	global $parser;

	list($chr, $start) = explode(" ", strtr($pos, ":", " "));
	
	list($stdout) = $parser->exec(get_path("samtools"), "depth -Q 1 -q 15 -r $chr:$start-$start $bam", true);
	$depth = trim(implode("", $stdout));
	if ($depth=="") return 0;
	
	return explode("\t", $depth)[2];
}

//returns a variant property (and optionally formats numbers)
function get_prop($var, $name, $digits = null)
{
	if (!isset($var[$name])) return "";
	if (!is_null($digits)) number_format($var[$name], $digits);
	return $var[$name];
}

//init
$ngsbits = get_path("ngs-bits");
$vcflib = get_path("vcflib");
$genome = genome_fasta($build);
$giab_bed = get_path("data_folder")."/dbs/GIAB/{$ref_sample}/high_conf_regions.bed";
if (!file_exists($giab_bed)) trigger_error("GiaB {$ref_sample} BED file missing: {$giab_bed}", E_USER_ERROR);
$giab_vcfgz = get_path("data_folder")."/dbs/GIAB/{$ref_sample}/high_conf_variants.vcf.gz";
if (!file_exists($giab_vcfgz)) trigger_error("GiaB {$ref_sample} VCF file missing: {$giab_vcfgz}", E_USER_ERROR);

//Target region base statistics
print "##Target region     : $roi\n";
$bases = get_bases($roi);
print "##Bases             : $bases\n";
//sort and merge $roi_hc after intersect - MH
$roi_hc = $parser->tempFile(".bed");
$pipeline = [];
$pipeline[] = ["{$ngsbits}BedIntersect", "-in $roi -in2 {$giab_bed}"];
$pipeline[] = ["{$ngsbits}BedSort", ""];
$pipeline[] = ["{$ngsbits}BedMerge", "-out $roi_hc"];
$parser->execPipeline($pipeline, "high-conf ROI");
$bases_hc = get_bases($roi_hc);
print "##High-conf bases   : $bases_hc (".number_format(100*$bases_hc/$bases, 2)."%)\n";
$roi_used = $roi_hc;
$bases_used = $bases_hc;
if ($min_dp>0)
{
	$roi_low_dp = $parser->tempFile(".bed");
	exec2("{$ngsbits}BedLowCoverage -bam {$bam} -in {$roi_hc} -cutoff {$min_dp} -out {$roi_low_dp}");
	$roi_high_dp = $parser->tempFile(".bed");
	exec2("{$ngsbits}BedSubtract -in {$roi_hc} -in2 {$roi_low_dp} -out {$roi_high_dp}");
	
	$bases_high_dp = get_bases($roi_high_dp);
	print "##High-depth bases  : $bases_high_dp (".number_format(100*$bases_high_dp/$bases, 2)."%)\n";
	
	$roi_used = $roi_high_dp;
	$bases_used = $bases_high_dp;
}
print "##Notice: Reference variants in the above region are evaluated!\n";
print "##\n";

//get reference variants in ROI
print "##Variant list      : $vcf\n";
$found = get_variants($vcf, $roi_used, false, $max_indel);
print "##Variants observed : ".count($found)."\n";
$expected = get_variants($giab_vcfgz, $roi_used, true, $max_indel);
print "##Variants expected : ".count($expected)."\n";

//find missing reference variants and variants with genotype mismatch
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
	else
	{
		$var_diff[$pos] = array("=", $var, $found[$pos]);
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

//sort differential variant list
function pos_gt($a, $b)
{
	list($c_a, $p_a) = explode(" ", strtr($a, array(":"=>" ")));
	list($c_b, $p_b) = explode(" ", strtr($b, array(":"=>" ")));
	
	if ($c_a==$c_b)
	{
		if ($p_a==$p_b) return 0;
		return ($p_a<$p_b) ? -1 : 1;
	}
	
	return ($c_a<$c_b) ? -1 : 1;
}
uksort($var_diff, "pos_gt");

//print TSV output
$sample = substr(basename($vcf), 0, 10);
print "#sample\tstatus\tpos\tvariant\tvariant_type\tref_GT\tobs_GT\tobs_DP\tobs_QUAL\tobs_MQM\tobs_AO\tobs_AF\tobs_QUAL_per_DP\n";
foreach($var_diff as $pos => list($status, $ref, $obs))
{
	if($status=="=" && !$matches) continue;
	
	$exchange = explode(" ", $pos)[1];
	$variant_type = strlen($exchange)==3 ? "SNV" : "INDEL";

	$ao = get_prop($obs, "AO");
	$dp = get_prop($obs, "DP");
	$gq = get_prop($obs, "GQ");
	$af = "n/a";
	if ($ao!="" && $dp!="" && $dp>0) $af = number_format($ao/$dp,2);
	
	$qual = get_prop($obs, "QUAL", 2);
	$qpd = "n/a";
	if ($qual!="" && $dp!="" && $dp>0) $qpd = number_format($qual/$dp, 2);
	
	print "$sample\t$status\t".strtr($pos, " ", "\t")."\t".$variant_type."\t".get_prop($ref, "GT")."\t".get_prop($obs, "GT")."\t".$dp."\t".$qual."\t".get_prop($obs, "MQM", 2)."\t".$ao."\t".$af."\t".$qpd."\n";
}
print "\n";


//print statistics
function stats_line($name, $options, $date, $expected, $var_diff, $var_type=null)
{
	global $bases_used;
	
	//count expected variants
	$c_expected = 0;
	foreach($expected as $pos => $var)
	{
		if (!is_null($var_type))
		{
			if (isset($var['TYPE']) && $var['TYPE']!=$var_type) continue;
		}
		
		++$c_expected;
	}
	
	//count FP/FN
	$fp = 0;
	$fn = 0;
	$gt_false = 0;
	foreach($var_diff as $pos => list($type, $ref, $obs))
	{
		//skip matches
		if ($type=="=") continue;
		
		if (!is_null($var_type))
		{
			if (isset($ref['TYPE']) && $ref['TYPE']!=$var_type) continue;
			if (isset($obs['TYPE']) && $obs['TYPE']!=$var_type) continue;
		}
		
		if ($type=="-")
		{
			$fn++;
		}
		if ($type=="g")
		{
			++$gt_false;
		}
		if ($type=="+")
		{	
			$fp++;
		}
	}
		
	$tp = $c_expected - $fn;
	$tn = $bases_used - $c_expected - $fp;
	$recall = number_format($tp/($tp+$fn),5);
	$precision = number_format($tp/($tp+$fp),5);
	$spec = number_format($tn/($tn+$fp),5);
	$geno_acc = number_format(($tp-$gt_false)/$tp,5);
	
	if (!is_null($var_type)) $name .= " ".$var_type;
	return implode("\t", array($name, $options, $date, $c_expected, $tp, $tn, $fp, $fn, $gt_false, $recall, $precision, $spec, $geno_acc));
}

if ($name=="") $name = basename($vcf, ".vcf.gz");
$options = array();
if ($min_dp>0) $options[] = "min_dp={$min_dp}";
if ($max_indel>0) $options[] = "max_indel={$max_indel}";
$options = implode(" ", $options);
$date = strtr(date("Y-m-d H:i:s", filemtime($vcf)), "T", " ");
$output = array();
$output[] = "#name\toptions\tdate\texpected variants\tTP\tTN\tFP\tFN\tgenotyping_errors\trecall/sensitivity\tprecision/ppv\tspecificity\tgenotyping_accuracy";
if (file_exists($stats))
{
	foreach(file($stats) as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		$output[] = $line;
	}
}
$output[] = stats_line($name, $options, $date, $expected, $var_diff);
$output[] = stats_line($name, $options, $date, $expected, $var_diff, "SNVS");
$output[] = stats_line($name, $options, $date, $expected, $var_diff, "INDELS");
file_put_contents($stats, implode("\n", $output)."\n");

?>
