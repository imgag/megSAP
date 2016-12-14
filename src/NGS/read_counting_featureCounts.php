<?php

/**
	@page read_counting_featureCounts
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("read_counting_featureCounts", "Perform read counting for aligned reads using featureCount contained in the subread package.");
$parser->addInfile("in",  "Input file in bam format.", false, "bam");
$parser->addString("prefix", "String to add to output files. Might include sub-directories.", false, NULL);

//optional parameters
$parser->addString("gtfFile", "GTF File containing feature annotations used for read counting.", true, get_path("data_folder")."genomes/gtf/ucsc_refseq_hg19.gtf");
$parser->addString("featureType", "Feature type used for mapping reads to features.", true, "exon");
$parser->addString("gtfAttribute", "GTF attribute used as feature ID.", true, "gene_id");
$parser->addFlag("stranded", "Specify whether a stranded protocol was used during library preparation. Default is non-stranded.");
$parser->addInt("minAQual", "Minimal alignment quality. Skip all reads with alignment quality lower than the given value. For paired-end reads at least one end should satisfy this criteria.", true, 3);
$parser->addFlag("paired", "The data is paired-end. Only properly paired reads are used.");
$parser->addFlag("includeMultiOverlap", "Count also reads that overlap more than one feature.");
$parser->addInt("threads", "Number of threads used for read counting", true, "4");

extract($parser->parse($argv));

//extracting sub-directories and generating folder structure
$out = create_path($prefix);
$sampleName = $out[0];
$outdir = $out[1];

$parser->log("read_counting_featureCounts output directory=$outdir");

//build command
$arguments = array();

$arguments[] = "$in";
$arguments[] = "-a $gtfFile";
$arguments[] = "-t $featureType";
$arguments[] = "-g $gtfAttribute";
$arguments[] = "-Q $minAQual";
$arguments[] = "-T $threads";

if($paired) {
	$arguments[] = "-p -B";
}

if($stranded) {
	$arguments[] = "-s 1";
} else {
	$arguments[] = "-s 0";
}

if($includeMultiOverlap) {
	$arguments[] = "-M";
}

$arguments[] = "-o ${outdir}${sampleName}_counts.tsv"; //output file

//execute command
$parser->log("Starting read counting");
$parser->exec(get_path("feature_counts"), implode(" ", $arguments), true);

// cleanup
$parser->exec("rm ", "${outdir}${sampleName}_counts.tsv.summary", false);
?>