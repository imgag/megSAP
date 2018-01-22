<?php

/**
 * @page vcf_somatic_donor
 */


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("vcf_somatic_donor", "Add donor germline variant information to somatic VCF files.");
$parser->addInfile("in_somatic", "Input somatic variant file (VCF).", false);
$parser->addInfileArray("in_donor", "Donor sample BAM file(s).", false);
$parser->addOutfile("out_vcf", "Output somatic variant file (VCF).", false);
extract($parser->parse($argv));

$tmp_vaf_in = $in_somatic;
$tmp_vaf_out = $parser->tempFile(".vcf");

foreach ($in_donor as $donor_bam)
{
	$args = [
		"-in", $tmp_vaf_in,
		"-bam", $donor_bam,
		"-out", $tmp_vaf_out,
		"-depth",
		"-name", "donor_" .basename($donor_bam, ".bam")
	];

	// add frequency and depth to somatic VCF file
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", implode(" ", $args), true);

	$tmp_vaf_in = $tmp_vaf_out;
	$tmp_vaf_out = $parser->tempFile(".vcf");
}

$parser->copyFile($tmp_vaf_in, $out_vcf);