<?php

/**
	@page spladder
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("spladder", "Run SplAdder on RNA-seq alignments for alternative splicing analysis.");

//mandatory arguments
$parser->addInfileArray("bams", "Input RNA-seq alignment(s) in BAM format", false);
$parser->addString("outdir", "Output directory", false);

//optional arguments
$parser->addString("gtf", "Annotation in GTF/GFF3 format", true, get_path("data_folder") . "/dbs/gene_annotations/GRCh37.gtf");
$parser->addInt("threads", "Number of threads", true, 2);

//extract arguments
extract($parser->parse($argv));

//run command
$arguments = array(
	"--bams=" . implode(",", $bams),
	"--outdir={$outdir}",
	"--annotation={$gtf}",
	#"--pyproc=y",
	"--parallel={$threads}",
	"--verbose=y"
);
$parser->exec(get_path("spladder"), implode(" ", $arguments), true);
