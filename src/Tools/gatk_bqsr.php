<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//QBSR docu
//https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
//https://gatk.broadinstitute.org/hc/en-us/articles/360035890531
//https://gatk.broadinstitute.org/hc/en-us/articles/5358896138011-BaseRecalibrator
//https://gatk.broadinstitute.org/hc/en-us/articles/5358826654875-ApplyBQSR

// parse command line arguments
$parser = new ToolBase("gatk_bqsr", "Calculates GATK BQSR model and applies it.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addOutfile("out",  "Output BAM file.", false);
//optional
$parser->addInt("threads",  "Number of threads to use.", true, 4);
$parser->addString("build", "Reference genome build.", true, "GRCh38");
extract($parser->parse($argv));

//init
$gatk = get_path("gatk");
$ref = genome_fasta($build);
$known_sites = "/mnt/storage2/GRCh38/share/data/dbs/gnomAD/gnomAD_genome_v3.1.2_GRCh38.vcf.gz";

//perform BQSR model generation
$bqsr_model = substr($out,0, -4)."_bqsr.recal_data_table";
$parser->exec($gatk, "BaseRecalibrator -I {$in} -R {$ref} --known-sites {$known_sites} -O {$bqsr_model}");

//apply model
$parser->exec($gatk, "ApplyBQSR -I {$in} -R {$ref} --bqsr-recal-file {$bqsr_model} -O {$out}");

?>