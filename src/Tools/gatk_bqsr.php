<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//BQSR docu
//https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
//https://gatk.broadinstitute.org/hc/en-us/articles/360035890531
//https://gatk.broadinstitute.org/hc/en-us/articles/5358896138011-BaseRecalibrator
//https://gatk.broadinstitute.org/hc/en-us/articles/5358826654875-ApplyBQSR

// parse command line arguments
$parser = new ToolBase("gatk_bqsr", "Calculates GATK BQSR model and applies it.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addOutfile("out",  "Output BAM file.", false);
//optional
$parser->addString("build", "Reference genome build.", true, "GRCh38");
extract($parser->parse($argv));

//init
/* $gatk = get_path("gatk"); */

$out_files = array();
$ref = genome_fasta($build);
$known_sites = "/mnt/storage2/megSAP/data/dbs/gnomAD/gnomAD_genome_v4.1_GRCh38.vcf.gz";

//input files for gatk BaseRecalibrator
$in_files = array();
$in_files[] = $in;
$in_files[] = $ref;
$in_files[] = $known_sites;

//perform BQSR model generation
$bqsr_model = substr($out,0, -4)."_bqsr.recal_data_table";

//ouput files for gatk BaseRecalibrator
$out_files[] = $bqsr_model;
$parser->execSingularity("gatk", get_path("container_gatk"), "gatk", "BaseRecalibrator -I {$in} -R {$ref} --known-sites {$known_sites} -O {$bqsr_model}", $in_files, $out_files);
/* $parser->exec($gatk, "BaseRecalibrator -I {$in} -R {$ref} --known-sites {$known_sites} -O {$bqsr_model}"); */

//apply model
//input files for gatk ApplyBQSR
$in_files = [
    $in,
    $ref,
    $bqsr_model
];
//output file for gatk ApplyBQSR
$out_files = [$out];
$parser->execSingularity("gatk", get_path("container_gatk"), "gatk", "ApplyBQSR -I {$in} -R {$ref} --bqsr-recal-file {$bqsr_model} -O {$out}", $in_files, $out_files);
/* $parser->exec($gatk, "ApplyBQSR -I {$in} -R {$ref} --bqsr-recal-file {$bqsr_model} -O {$out}"); */

?>