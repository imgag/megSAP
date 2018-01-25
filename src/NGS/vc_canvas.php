<?php

/**
	@page vc_canvas
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_canvas", "Somatic CNV calling with Canvas.");

// mandatory arguments
$parser->addInfile("tumor_bam", "Tumor sample BAM file.", false);
$parser->addInfile("normal_bam", "Normal sample BAM file.", false);
$parser->addString("out_dir", "Output directory.", false);

// optional arguments
$parser->addString("canvas_ref", "Canvas genome reference directory.",
	true, "/mnt/users/ahmattj1/tmp/canvas-grch37");

// extract arguments
extract($parser->parse($argv));

if (! is_dir($out_dir))
{
	mkdir($out_dir);
}

// create manifest file
$manifest_tmp = $parser->tempFile("manifest.txt");
$parser->exec("awk",
	"-v OFS=$'\t' '{ print $1\"_\"$2\"_\"$3,$1,$2 + 1,$3,0,0 }' </mnt/share/data/enrichment/ssHAEv6_2017_01_05.bed awk >${manifest_tmp}");

$baf_vcf = dirname($tumor_bam) . "/" . basename($tumor_bam, ".bam") . "_var.vcf.gz";

// run command
$arguments = [
	"Tumor-normal-enrichment",
	"-b", $tumor_bam,
	"--normal-bam", $normal_bam,
	"--reference", $canvas_ref . "kmer.fa",
	"--manifest", $manifest_tmp,
	"-g", $canvas_ref,
	"-n", basename($tumor_bam, ".bam"),
	"-f", $canvas_ref . "filter13.bed",
	"-o", $out_dir,
	"--b-allele-vcf", $baf_vcf,
	"--custom-parameters=CanvasBin,-m=TruncatedDynamicRange"
];
$parser->exec("/mnt/share/opt/Canvas-1.31.0.843+master_x64/Canvas",
	implode(" ", $arguments), true);