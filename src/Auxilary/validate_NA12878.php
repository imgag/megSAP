<?php
/** 
	@page validate_NA12878
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("validate_NA12878", "Validates the performance of a sequencing experiment on the GiaB reference sample NA12878.");
$parser->addInfile("vcf", "Variant list (VCF or VCF.GZ).", false);
$parser->addInfile("bam", "Mapped reads (BAM format).", false);
$parser->addInfile("roi", "Target region of the processing system (BED format).", false);
$parser->addOutfile("stats", "Append statistics to this file.", false);
//optional
$parser->addString("name", "Name used in the 'stats' output. If unset, the 'vcf' file base name is used.", true);
$parser->addString("sample_id", "Sample ID in VCF heanders (needed only if several samples are present in the VCF, e.g. in a trio analysis.", true, "");
$parser->addInt("min_dp", "If set, only regions in the 'roi' with at least the given depth are evaluated.", true, 0);
$parser->addInt("min_qual", "If set, only input variants with QUAL greater or equal to the given value are evaluated.", true, 5);
$parser->addInt("max_indel", "Maximum indel size (larger indels are ignored). Disabled if set to 0.", true, 50);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("ref_sample", "Reference sample to use for validation.", true, "NA12878");
$parser->addFlag("matches", "Do not only show variants that were missed (-), are novel (+) or with genotype mismatch (g), but also show matches (=) in output.");
$parser->addFlag("skip_depth_calculation", "Do not calculate depth of missed variants to speed up calculation.");
extract($parser->parse($argv));

//returns the variants of a VCF file in the given ROI
function get_variants($vcf_gz, $roi, $max_indel, $min_qual, $sample_id, &$skipped)
{
	global $parser;
	global $genome;
		
	//get variants
	$tmp = $parser->tempFile(".vcf");
	$pipeline = [];
	$pipeline[] = array(ends_with(strtolower($vcf_gz), ".vcf") ? "cat" : "zcat", $vcf_gz);
	$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfFilter", "-reg {$roi} -ref $genome", [$roi, $genome], [], true));
	$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfStreamSort", "-out {$tmp}", [], [], true));
	$parser->execPipeline($pipeline, "variant extraction");
	
	//put together output
	$output = array();	
	$file = file($tmp);
	$is_trio = false;
	foreach($file as $line)
	{
		$line = trim($line);

		if (starts_with(strtolower($line), "#chrom"))
		{
			$header = explode("\t", $line);
			if (count($header) > 10) 
			{
				$is_trio = true;
				$idx_child = array_search($sample_id, $header);
				if ($idx_child == false) trigger_error("Column for child ID '$sample_id' could not be found in $vcf_gz.", E_USER_ERROR);
			}
		}

		if ($line=="" || $line[0]=="#") continue;

		//get variant infos
		$columns = explode("\t", $line);
		if ($is_trio)
		{
			list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = array_slice($columns, 0, 9);
			$sample = $columns[$idx_child];
		}
		else list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = $columns;
		
		if (!starts_with($chr, "chr")) $chr = "chr".$chr;
		
		
		//multi-allelic > abort
		if (contains($alt, ",")) trigger_error("This tool cannot handle multi-allelic variants. Break and normalize variants using VcfBreakMulti, VcfLeftNormalize and VcfStreamSort before running this tool!", E_USER_ERROR);
		
		//skip low mappabilty variants (megSAP)
		if (contains($filter, "low_mappability"))
		{
			@$skipped["low_mappability"] += 1;
			continue;
		}
		
		//skip mosaic variants (megSAP)
		if (contains($filter, "mosaic"))
		{
			@$skipped["mosaic"] += 1;
			continue;
		}
		
		//skip low mappabilty variants (DRAGEN)
		//TODO
		
		//skip mosaic variants (DRAGEN)
		$info = explode(";",$info);
		foreach($info as $entry)
		{
			if ($entry=="MOSAIC")
			{
				@$skipped["mosaic"] += 1;
				continue(2);
			}
		}
		
		//compile output
		$var = array();
		$var["QUAL"] = $qual;
		foreach($info as $entry)
		{
			if (starts_with($entry, "MQM="))
			{
				$var['MQM'] = substr($entry, 4);
			}
		}
		
		$format = explode(":", $format);
		$sample = explode(":", $sample);
		
		//get GT from FORMAT
		if ($format[0]!="GT") trigger_error("First element of FORMAT is not GT in $line", E_USER_ERROR);
		$gt = strtr($sample[0], array("|"=>"/", "."=>"0"));
		if ($gt=="1/0") $gt="0/1";
		$var['GT'] = $gt;
		
		//skip wildtype variants (should not happen in normal pipeline, but can happen)
		if ($gt=="0/0")
		{
			@$skipped["genotype wildtype: 0/0"] += 1;
			continue;
		}
		if ($gt!="0/1" && $gt!="1/1")
		{
			@$skipped["genotype invalid: ".$sample[0]] += 1;
			continue;
		}
		
		//get rest from FORMAT
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
			if ($max_indel>0 && (strlen($ref)>$max_indel || strlen($alt)>$max_indel))
			{
				@$skipped["indel>{$max_indel}"] += 1;
				continue;
			}
		}
		else
		{
			$var["TYPE"] = "SNVS";
		}
		
		//filter by QUAL
		if ($min_qual>0 && $qual<$min_qual)
		{
			@$skipped["QUAL<{$min_qual}"] += 1;
			continue;
		}
		
		$tag = "{$chr}:{$pos} {$ref}>{$alt}";
		if (isset($output[$tag]))
		{
			@$skipped["variant found more than once"] += 1;
			continue;
		}
		$output[$tag] = $var;
	}
	return $output;
}

//returns the depth at a specified position in the BAM file
function get_depth($pos, $bam)
{
	global $parser;

	list($chr, $start) = explode(" ", strtr($pos, ":", " "));
	
	list($stdout) = $parser->execApptainer("samtools", "samtools depth", "-Q 1 -q 15 -r $chr:$start-$start $bam", [$bam]);
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
$genome = genome_fasta($build);
$giab_bed = get_path("data_folder")."/dbs/GIAB/{$ref_sample}/high_conf_regions.bed";
if (!file_exists($giab_bed)) trigger_error("GiaB {$ref_sample} BED file missing: {$giab_bed}", E_USER_ERROR);
$giab_vcfgz = get_path("data_folder")."/dbs/GIAB/{$ref_sample}/high_conf_variants_normalized.vcf.gz";
if (!file_exists($giab_vcfgz)) trigger_error("GiaB {$ref_sample} VCF file missing: {$giab_vcfgz}", E_USER_ERROR);

//Target region base statistics
print "##Target region     : $roi\n";
$bases = bed_size($roi);
print "##Bases             : $bases\n";
//sort and merge $roi_hc after intersect - MH
$roi_hc = $parser->tempFile(".bed");
$pipeline = [];
$pipeline[] = array("", $parser->execApptainer("ngs-bits", "BedIntersect", "-in $roi -in2 {$giab_bed}", [$roi, $giab_bed], [], true));
$pipeline[] = array("", $parser->execApptainer("ngs-bits", "BedSort", "", [], [], true));
$pipeline[] = array("", $parser->execApptainer("ngs-bits", "BedMerge", "-out $roi_hc", [], [], true));
$parser->execPipeline($pipeline, "high-conf ROI");
$bases_hc = bed_size($roi_hc);
print "##High-conf bases   : $bases_hc (".number_format(100*$bases_hc/$bases, 2)."%)\n";
$roi_used = $roi_hc;
$bases_used = $bases_hc;
if ($min_dp>0)
{
	$roi_low_dp = $parser->tempFile(".bed");
	$parser->execApptainer("ngs-bits", "BedLowCoverage", "-bam {$bam} -in {$roi_hc} -cutoff {$min_dp} -out {$roi_low_dp} -threads 4 -ref {$genome}", [$bam, $genome]);
	$roi_high_dp = $parser->tempFile(".bed");
	$parser->execApptainer("ngs-bits", "BedSubtract", "-in {$roi_hc} -in2 {$roi_low_dp} -out {$roi_high_dp}");
	
	$bases_high_dp = bed_size($roi_high_dp);
	print "##High-depth bases  : $bases_high_dp (".number_format(100*$bases_high_dp/$bases, 2)."%)\n";
	
	$roi_used = $roi_high_dp;
	$bases_used = $bases_high_dp;
}
print "##Notice: Reference variants in the above region are evaluated!\n";
print "##\n";

//get reference variants in ROI
print "##Variant list      : $vcf\n";
$skipped = [];
$found = get_variants($vcf, $roi_used, $max_indel, $min_qual, $sample_id, $skipped);
print "##Variants observed : ".count($found)."\n";
foreach($skipped as $reason => $count)
{
	if ($count==0) continue;
	print "##  Skipped '{$reason}' variants: {$count}\n";
}
$skipped = [];
$expected = get_variants($giab_vcfgz, $roi_used, $max_indel, 0, "", $skipped);
print "##Variants expected : ".count($expected)."\n";
foreach($skipped as $reason => $count)
{
	if ($count==0) continue;
	print "##  Skipped '{$reason}' variants: {$count}\n";
}

//find missing reference variants and variants with genotype mismatch
$var_diff = array();
foreach($expected as $pos => $var)
{
	if (!isset($found[$pos]))
	{
		$exp = array();
		if($skip_depth_calculation)
		{
			$exp["DP"] = "N/A";
		}
		else
		{
			$exp["DP"] = get_depth($pos, $bam);
		}
		
		
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
	
	print "$name\t$status\t".strtr($pos, " ", "\t")."\t".$variant_type."\t".get_prop($ref, "GT")."\t".get_prop($obs, "GT")."\t".$dp."\t".$qual."\t".get_prop($obs, "MQM", 2)."\t".$ao."\t".$af."\t".$qpd."\n";
}
print "\n";


//print statistics
function stats($var_type, $expected, $var_diff)
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
	$recall = number_format($tp/($tp+$fn),4);
	$precision = number_format($tp/($tp+$fp),4);
	$geno_acc = number_format(($tp-$gt_false)/$tp,4);
	$f1 = number_format(2*$tp/(2*$tp+$fp+$fn),4);
	
	return [$c_expected, $recall, $precision, $geno_acc, $f1];
}

if ($name=="") $name = basename($vcf, ".vcf.gz");
$options = array();
if ($min_dp>0) $options[] = "min_dp={$min_dp}";
if ($max_indel>0) $options[] = "max_indel={$max_indel}";
if ($min_qual>0) $options[] = "min_qual={$min_qual}";
$options = implode(" ", $options);
$date = strtr(date("Y-m-d H:i:s", filemtime($vcf)), "T", " ");
$output = array();
$output[] = "#name\toptions\tdate\taverage_depth\t20x percentage\texpected_snvs\texpected_indels\tsnv_sensitivity\tsnv_ppv\tsnv_f1\tsnv_genotyping_accuracy\tindel_sensitivity\tindel_ppv\tindel_f1\tindel_genotyping_accuracy\tall_sensitivity\tall_ppv\tall_f1\tall_genotyping_accuracy";
if (file_exists($stats))
{
	foreach(file($stats) as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		$output[] = $line;
	}
}

//get QC metrics
$qcml = dirname($bam)."/".basename2($bam)."_stats_map.qcML";
$avg_depth = get_qcml_value($qcml, "QC:2000025", "n/a");
$cov20x = get_qcml_value($qcml, "QC:2000027", "n/a");

list($snv_exp, $snv_sens, $snv_ppv, $snv_geno, $snv_f1) = stats("SNVS", $expected, $var_diff);
list($indel_exp, $indel_sens, $indel_ppv, $indel_geno, $indel_f1) = stats("INDELS", $expected, $var_diff);
list($all_exp, $all_sens, $all_ppv, $all_geno, $all_f1) = stats(null, $expected, $var_diff);
$output[] = implode("\t", [$name, $options, $date, $avg_depth, $cov20x, $snv_exp, $indel_exp, $snv_sens, $snv_ppv, $snv_f1, $snv_geno, $indel_sens, $indel_ppv, $indel_f1, $indel_geno, $all_sens, $all_ppv, $all_f1, $all_geno])."\n";
file_put_contents($stats, implode("\n", $output)."\n");
?>
