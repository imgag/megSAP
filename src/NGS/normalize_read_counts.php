<?php
/**
 * @page normalize_read_counts
 */

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("normalize_read_counts", "\$Rev: 2$", "Normalizes read counts produced by HTSeq using the topas toolkit");
$parser->addInfile("in", "Input read count file in tsv format", false, true);
$parser->addOutfile("out", "Output file for normalized read counts", false);

//optional parameters
$parser->addString("gtf", "Reference annotation used for read counting in GTF format", true, get_path("data_folder")."genomes/gtf/ucsc_refseq_hg19.gtf");
$parser->addString("method", "The normalization method. Options are: 'fpkm' or 'cpm'", true, "fpkm");
$parser->addString("feature", "The GTF feature used for read counting", true, "exon");
$parser->addString("idattr", "The ID attribute used to summarize feature types during read counting", true, "gene_id");
$parser->addFlag("header", "Indicate whether the input file has a header or not");

extract($parser->parse($argv));

//extracting sub-directories and generating folder structure
create_path($out);

$arguments = array();
$arguments[] = "NormExprTable";
$arguments[] = "-i $in";
$arguments[] = "-o $out";
$arguments[] = "-gtf $gtf";

if($method == "fpkm") {
	$arguments[] = "-rpkm";
} else {
	$arguments[] = "-cpm";
}

$arguments[] = "-type $feature";
$arguments[] = "-idattr $idattr";

if($header) {
	$arguments[] = "-header true";
} else {
	$arguments[] = "-header false";
}

$parser->log("Starting read count normalization");
$parser->exec(get_path("topas"), implode(" ", $arguments), true);
?>