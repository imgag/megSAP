<?php
/** 
	@page validate_HCC1395-HCC1395BL
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("validate_HCC1395-HCC1395BL", "Validates the somatic variant calling performance of a specified caller on the tumor-normal pair HCC1395-HCC1395BL.");
$parser->addInfile("vcf", "Called somatic variant list to validate on (VCF or VCF.GZ).", false);
$parser->addInfile("roi", "Target region of the processing system (BED format).", false);
$parser->addInfile("bam", "Mapped tumor reads (BAM format).", false);
$parser->addOutfile("stats", "Append statistics to this file.", false);
$parser->addEnum("caller", "The caller used to call somatic variants", false, ["strelka2", "dragen", "deepsomatic", "varscan2"]);
//optional
$parser->addFlag("ignore_filters", "Ignores filter column entries.");
$parser->addFlag("with_low_evs", "Ignores 'LowEVS' filter column entry for strelka2 callings");
$parser->addFlag("matches", "Do not only show variants that were missed (-) or novel (+), but also show matches (=) in output.");
$parser->addFlag("skip_depth_calculation", "Do not calculate depth of missed variants to speed up calculation.");
$parser->addFlag("soft", "Set to enable softer benchmarking including variants flagged as MedConf in truth set");
$parser->addInt("min_dp", "If set, only regions in the 'roi' with at least the given depth are evaluated.", true, 0);
$parser->addInt("max_indel", "Maximum indel size (larger indels are ignored). Disabled if set to 0.", true, 50);
$parser->addString("ref_sample", "Reference sample to use for validation.", true, "HCC1395");
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("name", "Name used in the 'stats' output. If unset, the 'vcf' file base name is used.", true);

extract($parser->parse($argv));

function get_variants_som($vcf_gz, $roi, $max_indel, $caller, &$skipped)
{
	global $parser;
	global $genome;
	global $with_low_evs;
	global $ignore_filters;
	global $soft;

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
	foreach($file as $line)
	{
		$line = trim($line);

		if ($line=="" || $line[0]=="#") continue;

		if ($caller == "ref") list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = explode("\t", $line."\t");
		else if ($caller == "deepsomatic" || $caller == "varscan2") list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample_tumor) = explode("\t", $line."\t");
		else list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample_normal, $sample_tumor) = explode("\t", $line."\t");

		if (!is_numeric(chr_trim($chr))) 
		{
			@$skipped["non-autosomes"] += 1;
			continue; //skip non-autosomes
		}

		//fix chr
		if (!starts_with($chr, "chr")) $chr = "chr".$chr;

		//multi-allelic > abort
		if (contains($alt, ",")) trigger_error("This tool cannot handle multi-allelic variants. Break and normalize variants using VcfBreakMulti, VcfLeftNormalize and VcfStreamSort before running this tool!", E_USER_ERROR);

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

		if (strlen($ref)>1 || strlen($alt)>1)
		{
			$var["TYPE"] = "INDELS";
			if ($max_indel>0 && (strlen($ref)>$max_indel || strlen($alt)>$max_indel))
			{
				@$skipped["indel>{$max_indel}"] += 1;
				continue;
			}
		} 
		else $var["TYPE"] = "SNVS";

		if ($caller == "strelka2")
		{
			list($tumor_dp, $tumor_af) = ($var["TYPE"] == "INDELS") ? vcf_strelka_indel($format, $sample_tumor) : vcf_strelka_snv($format, $sample_tumor, $alt);
			
			//clear filter column
			$filter_replace = ["PASS"=>"", "freq-tum"=>"", ";"=>""];
			if ($with_low_evs) $filter_replace ["LowEVS"] = "";
		}
		else if ($caller == "dragen")
		{
			list($tumor_dp, $tumor_af) = vcf_dragen_var($format, $sample_tumor);
			
			//clear filter column
			$filter_replace = ["PASS"=>"", "."=>"", "freq-tum"=>"", ";"=>""];
		}
		else if ($caller == "deepsomatic")
		{
			list($tumor_dp, $tumor_af) = vcf_deepvariant($format, $sample_tumor);

			//clear filter column
			$filter_replace = ["PASS"=>"", "freq-tum"=>"", ";"=>""];
		}
		else if ($caller == "varscan2")
		{
			list($tumor_dp, $tumor_af) = vcf_varscan2($format, $sample_tumor);
			//clear filter column
			$filter_replace = ["PASS"=>"", "."=>"", "freq-tum"=>"", ";"=>""];
		}
		else if ($caller = "ref")
		{
			//clear filter column
			$filter_replace = ["PASS"=>"", "HighConf"=>"", ";"=>""];
			$filter = trim(strtr($filter, $filter_replace));
			if ($soft) $filter = "";
			
			if ($filter != "") 
			{
				@$skipped["MedConf"] += 1;
				continue;
			}
		}
		else
		{
			trigger_error("Unknown 'caller' given in load_vcf(): $caller", E_USER_ERROR);
		}
		
		if ($caller != "ref")
		{
			$filter = trim(strtr($filter, $filter_replace));
			if ($ignore_filters) $filter = "";

			if ($filter != "") 
			{
				@$skipped["filtered"] += 1;
				continue;
			}
			
			$var["DP"] = $tumor_dp;
			$var["AF"] = $tumor_af;
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

//print statistics
function stats_som($var_type, $expected, $var_diff)
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
		if ($type=="+")
		{	
			$fp++;
		}
	}
		
	$tp = $c_expected - $fn;
	$tn = $bases_used - $c_expected - $fp;
	$recall = number_format($tp/($tp+$fn),4);
	$precision = number_format($tp/($tp+$fp),4);
	$f1 = number_format(2*$tp/(2*$tp+$fp+$fn),4);
	
	return [$c_expected, $recall, $precision, $f1];
}

//init
if ($name=="") $name = basename($vcf, ".vcf.gz");
$genome = genome_fasta($build);

$hc_bed = get_path("data_folder")."/dbs/GIAB/{$ref_sample}/high_conf_regions.bed";
if (!file_exists($hc_bed)) trigger_error("High-confidence BED for {$ref_sample} missing: {$hc_bed}", E_USER_ERROR);
$truth_vcf = get_path("data_folder")."/dbs/GIAB/{$ref_sample}/high_conf_variants_normalized.vcf.gz";
if (!file_exists($truth_vcf)) trigger_error("{$ref_sample} VCF file missing: {$truth_vcf}", E_USER_ERROR);

//Target region base statistics
print "##Target region     : $roi\n";
$bases = bed_size($roi);
print "##Bases             : $bases\n";
//sort and merge $roi_hc after intersect - MH
$roi_hc = $parser->tempFile(".bed");
$pipeline = [];
$pipeline[] = array("", $parser->execApptainer("ngs-bits", "BedIntersect", "-in $roi -in2 {$hc_bed}", [$roi, $hc_bed], [], true));
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
$found = get_variants_som($vcf, $roi_used, $max_indel, $caller, $skipped);
print "##Variants observed : ".count($found)."\n";
foreach($skipped as $reason => $count)
{
	if ($count==0) continue;
	print "##  Skipped '{$reason}' variants: {$count}\n";
}
$skipped = [];
$expected = get_variants_som($truth_vcf, $roi_used, $max_indel, "ref", $skipped);
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
		if($skip_depth_calculation) $exp["DP"] = "N/A";
		else $exp["DP"] = get_depth($pos, $bam);		
		
		$var_diff[$pos] = array("-", $var, $exp);
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

//print TSV output
print "#sample\tstatus\tpos\tvariant\tvariant_type\tobs_DP\tobs_QUAL\tobs_MQM\tobs_AF\tobs_QUAL_per_DP\n";
foreach($var_diff as $pos => list($status, $ref, $obs))
{
	if($status=="=" && !$matches) continue;
	
	$exchange = explode(" ", $pos)[1];

	$variant_type = strlen($exchange)==3 ? "SNV" : "INDEL";
	$af = get_prop($obs, "AF");
	$dp = get_prop($obs, "DP");
	
	$qual = get_prop($obs, "QUAL");
	$qpd = "n/a";
	if (is_numeric($qual) && is_numeric($dp) && $dp > 0) $qpd = number_format(floatval($qual) / floatval($dp), 2);
	
	print "$name\t$status\t".strtr($pos, " ", "\t")."\t".$variant_type."\t".$dp."\t".$qual."\t".get_prop($obs, "MQM", 2)."\t".$af."\t".$qpd."\n";
}
print "\n";

//print statistics
$options = array();
if ($min_dp>0) $options[] = "min_dp={$min_dp}";
if ($max_indel>0) $options[] = "max_indel={$max_indel}";
if ($soft) $options[] = "soft_benchmark";
if ($ignore_filters) $options[] = "ignore_filters";
$options = implode(" ", $options);

$date = strtr(date("Y-m-d H:i:s", filemtime($vcf)), "T", " ");
$output = array();
$output[] = "#name\toptions\tdate\taverage_depth\t20x percentage\texpected_snvs\texpected_indels\tsnv_sensitivity\tsnv_ppv\tsnv_f1\tindel_sensitivity\tindel_ppv\tindel_f1\tall_sensitivity\tall_ppv\tall_f1";
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
$avg_depth = "n/a";
if (file_exists($qcml))
{
	list($stdout) = exec2("grep 'QC:2000025' $qcml"); //Average sequencing depth in target region
	foreach(explode(" ", $stdout[0]) as $part)
	{
		if (!contains($part, "value")) continue;
		$avg_depth = explode("\"", $part)[1];
	}
}
$cov20x = "n/a";
if (file_exists($qcml))
{
	list($stdout) = exec2("grep 'QC:2000027' $qcml"); //Average sequencing depth in target region
	foreach(explode(" ", $stdout[0]) as $part)
	{
		if (!contains($part, "value")) continue;
		$cov20x = explode("\"", $part)[1];
	}
}

list($snv_exp, $snv_sens, $snv_ppv, $snv_f1) = stats_som("SNVS", $expected, $var_diff);
list($indel_exp, $indel_sens, $indel_ppv, $indel_f1) = stats_som("INDELS", $expected, $var_diff);
list($all_exp, $all_sens, $all_ppv, $all_f1) = stats_som(null, $expected, $var_diff);
$output[] = implode("\t", [$name, $options, $date, $avg_depth, $cov20x, $snv_exp, $indel_exp, $snv_sens, $snv_ppv, $snv_f1, $indel_sens, $indel_ppv, $indel_f1, $all_sens, $all_ppv, $all_f1])."\n";
file_put_contents($stats, implode("\n", $output)."\n");
?>
