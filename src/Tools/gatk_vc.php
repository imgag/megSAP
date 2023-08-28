<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//HaplotypeCaller docu
//https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode
//https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
//https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
//https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle > https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

// parse command line arguments
$parser = new ToolBase("gatk_vc", "Small variant calling with GATK HaplotypeCaller.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addInfile("roi", "Input BED file with target region for calling.", false);
$parser->addOutfile("out",  "Output VCF.GZ file (indexed).", false);
//optional
$parser->addFlag("pcrfree", "Do not use PRC model for InDels.");
$parser->addFlag("gvcf", "Produce gVCF instead of normal VCF (also skips post-processing).");
$parser->addString("build", "Reference genome build.", true, "GRCh38");
$parser->addInt("threads", "Threads to use (attention: scaling is poor).", true, 4);
extract($parser->parse($argv));

//init
$gatk = get_path("gatk");
$ref = genome_fasta($build);

//perform variant calling
$args = [];
$args[] = "-L ".$roi;
$args[] = "-G StandardAnnotation";
$args[] = "-G StandardHCAnnotation";
$args[] = "--native-pair-hmm-threads ".$threads;
if ($pcrfree)
{
	$args[] = "--pcr-indel-model NONE";	
}
if ($gvcf)
{
	$args[] = "-ERC GVCF";
	$args[] = "-G AS_StandardAnnotation";
}

$tmp = $gvcf ? $out : $parser->tempFile(".vcf.gz");
$parser->exec($gatk, "HaplotypeCaller -R {$ref} -I {$in} -O {$tmp} ".implode(" ", $args));

if (!$gvcf)
{
	//perform postprocessing
	$pipeline = [];
	$pipeline[] = array("zcat", $tmp);
	$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");
	$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize", "-stream -ref {$ref}");
	$pipeline[] = array(get_path("ngs-bits")."VcfStreamSort", "");
	$pipeline[] = array("bgzip", "-c > {$out}", false);
	$parser->execPipeline($pipeline, "post processing");

	//index VCF
	$parser->exec("tabix", "-p vcf {$out}", false);
}

?>