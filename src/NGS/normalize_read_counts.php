<?php
/**
 * @page normalize_read_counts
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("normalize_read_counts", "Normalizes read counts produced by HTSeq using the topas toolkit");
$parser->addInfile("in", "Input read count file in tsv format", false, true);
$parser->addOutfile("out", "Output file for normalized read counts", false);
//optional parameters
$parser->addString("gtf", "Reference annotation used for read counting in GTF format", true, get_path("data_folder")."genomes/gtf/ucsc_refseq_hg19.gtf");
$parser->addEnum("method", "The normalization method.", true, array("rpkm", "cpm"), "rpkm");
$parser->addString("feature", "The GTF feature used for read counting", true, "exon");
$parser->addString("idattr", "The ID attribute used to summarize feature types during read counting", true, "gene_id");
$parser->addFlag("header", "Indicate whether the input file has a header or not");
extract($parser->parse($argv));

$arguments = array();
$arguments[] = "NormExprTable";
$arguments[] = "-i $in";
$arguments[] = "-o $out";
$arguments[] = "-gtf $gtf";
$arguments[] = "-".$method;
$arguments[] = "-type $feature";
$arguments[] = "-idattr $idattr";
$arguments[] = "-header ".($header ? "true" : "false");
$parser->exec(get_path("topas"), implode(" ", $arguments), true);
?>