<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_cfdna", "cfDNA mutation caller.");
$parser->addInfile("bam", "Input BAM file, with deduplicated alignments and DP tags.", false);
$parser->addInfile("target", "Target region BED file.", false);
$parser->addOutfile("vcf", "Variant call output as VCF file.", true);
$parser->addOutfile("tsv", "Variant call output as TSV file.", true);
$parser->addOutfile("mrd", "MRD probability output file.", true);
$parser->addInfile("model", "Error model parameters.", true);
$parser->addInfile("model_region", "Target region to use for parameter estimation.", true);

extract($parser->parse($argv));

?>
