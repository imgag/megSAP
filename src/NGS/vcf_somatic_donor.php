<?php

/**
 * @page vcf_intersect
 */


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("vcf_somatic_donor", "Mark potential donor germline variants in somatic VCF files.");

$parser->addInfile("in_somatic", "Input somatic variant file (VCF).", false);
$parser->addInfile("in_donor", "Donor sample BAM file.", false);

$parser->addOutfile("out_vcf", "Output somatic variant file (VCF).", false);

//optional
$parser->addFloat("min_af", "Minimum variant allele frequency in donor germline analysis.", true, 0.05);
$parser->addFloat("min_depth", "Minimum variant depth donor germline analysis.", true, 20);

extract($parser->parse($argv));

// 1. VAF -in somatic.vcf -out tmp.vcf -depth -name donor
// 2. VCF filter with somatic_donor filter

$tmp_vaf_out = $parser->tempFile(".vcf");

$args = [
	"-in", $in_somatic,
	"-bam", $in_donor,
	"-out", $tmp_vaf_out,
	"-depth",
	"-name", "donor"
];

// add frequency and depth to somatic VCF file
$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", implode(" ", $args), true);

// add somatic-donor filter
$parser->execTool("NGS/filter_vcf.php", "-in {$tmp_vaf_out} -out {$out_vcf} -type somatic-donor -keep");