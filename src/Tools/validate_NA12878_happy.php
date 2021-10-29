<?php
/** 
	@page validate_NA12878_happy
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("validate_NA12878_happy", "Validates the performance of a sequencing experiment on the GiaB reference sample NA12878 using hap.py.");
$parser->addInfile("vcf", "Input variant list of sequencing experiment (VCF.GZ format). Note: Multi-allelic variants must be split and indels must be left-aligned.", false);
$parser->addInfile("roi", "Target region to evaluate (BED format).", false);
$parser->addOutfile("stats", "Append statistics to this file.", false);
//optional
$parser->addString("name", "Name used in the 'stats' output. If unset, the 'vcf' file base name is used.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("ref_sample", "Reference sample to use for validation.", true, "NA12878");
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

//returns the base count of a BED file
function get_variants($vcf)
{
	list($stdout) = exec2("cat $vcf | egrep -v \"^#\" | wc -l");
	return trim($stdout[0]);
}

//print statistics
function stats_line($name, $options, $date, $bases_hc, $query_total, $fp, $fn, $fp_gt)
{
	$tp = $query_total - $fn;
	$tn = $bases_hc - $query_total - $fp;
	$recall = number_format($tp/($tp+$fn),5);
	$precision = number_format($tp/($tp+$fp),5);
	$spec = number_format($tn/($tn+$fp),5);
	$geno_acc = number_format(($tp-$fp_gt)/$tp,5);
	
	return implode("\t", array($name, $options, $date, $query_total, $tp, $tn, $fp, $fn, $fp_gt, $recall, $precision, $spec, $geno_acc));
}

//init
$tmp_folder = $parser->tempFolder("validate_NA12878_happy");
$happy = get_path("happy");
$ngsbits = get_path("ngs-bits");
$genome = genome_fasta($build);
$giab_bed = get_path("data_folder")."/dbs/GIAB/{$ref_sample}/high_conf_regions.bed";
if (!file_exists($giab_bed)) trigger_error("GiaB {$ref_sample} BED file missing: {$giab_bed}", E_USER_ERROR);
$giab_vcfgz = get_path("data_folder")."/dbs/GIAB/{$ref_sample}/high_conf_variants.vcf.gz";
if (!file_exists($giab_vcfgz)) trigger_error("GiaB {$ref_sample} VCF file missing: {$giab_vcfgz}", E_USER_ERROR);

//target region handling
print "##Target region  : $roi\n";
$bases = get_bases($roi);
print "##Bases          : $bases\n";
//sort and merge $roi_hc after intersect - MH
$roi_hc = $tmp_folder."/roi_hc.bed";
$pipeline = [];
$pipeline[] = ["{$ngsbits}BedIntersect", "-in $roi -in2 {$giab_bed}"];
$pipeline[] = ["{$ngsbits}BedSort", ""];
$pipeline[] = ["{$ngsbits}BedMerge", "-out $roi_hc"];
$parser->execPipeline($pipeline, "high-conf ROI");
$bases_hc = get_bases($roi_hc);
print "##High-conf region: $roi_hc\n";
print "##High-conf bases: $bases_hc (".number_format(100*$bases_hc/$bases, 2)."%)\n";
print "##Notice: Reference variants in the above region are evaluated!\n";

//perform comparison
list($stdout, $stderr) = $parser->exec($happy, "{$giab_vcfgz} {$vcf} -T {$roi_hc} -r {$genome} -o {$tmp_folder}/output", true);

//append output to statistics file
if ($name=="") $name = basename($vcf, ".vcf.gz");
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

$all_truth_total = 0;
$all_fp = 0;
$all_fn = 0;
$all_fp_gt = 0;
$out_snv = "";
$out_indel = "";
foreach(file("{$tmp_folder}/output.summary.csv") as $line)
{
	if (starts_with($line, "Type")) continue;
	list($type, $filter, $truth_total, $tp, $fn, $query_total, $fp, $unk, $fp_gt, $fp_al) = explode(",", $line); 
	if($filter=="PASS") continue;
	
	$all_truth_total += $truth_total;
	$all_fp += $fp;
	$all_fn += $fn;
	$all_fp_gt += $fp_gt;
	if ($type=="SNP")
	{
		$out_snv = stats_line($name." SNVS", "", $date, $bases_hc, $truth_total, $fp, $fn, $fp_gt);
	}
	if ($type=="INDEL")
	{
		$out_indel = stats_line($name." INDELS", "", $date, $bases_hc, $truth_total, $fp, $fn, $fp_gt);
	}
}
$output[] = stats_line($name, "", $date, $bases_hc, $all_truth_total, $all_fp, $all_fn, $all_fp_gt);
$output[] = $out_snv;
$output[] = $out_indel;
file_put_contents($stats, implode("\n", $output)."\n");

?>
