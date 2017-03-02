<?php

/**
	@page read_counting_featureCounts
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("read_counting_featureCounts", "Perform read counting for aligned reads using featureCount contained in the subread package.");
$parser->addInfile("in",  "BAM input file.", false, "bam");
$parser->addOutfile("out", "Raw count output TSV file.", false);

//optional parameters
$parser->addString("gtfFile", "GTF File containing feature annotations used for read counting.", true, get_path("data_folder")."/dbs/UCSC/refGene.gtf");
$parser->addString("featureType", "Feature type used for mapping reads to features.", true, "exon");
$parser->addString("gtfAttribute", "GTF attribute used as feature ID.", true, "gene_id");
$parser->addFlag("stranded", "Specify whether a stranded protocol was used during library preparation. Default is non-stranded.");
$parser->addInt("minAQual", "Minimal alignment quality. Skip all reads with alignment quality lower than the given value. For paired-end reads at least one end should satisfy this criteria.", true, 3);
$parser->addFlag("paired", "The data is paired-end. Only properly paired reads are used.");
$parser->addFlag("includeMultiOverlap", "Count reads multiple times if they overlap more than one feature.");
$parser->addInt("threads", "Number of threads used for read counting", true, "4");
extract($parser->parse($argv));

//build arguments array
$args = array();
$args[] = "-a $gtfFile";
$args[] = "-t $featureType";
$args[] = "-g $gtfAttribute";
$args[] = "-Q $minAQual";
$args[] = "-T $threads";
$args[] = "-s ".($stranded ? "1" : "0");
if($paired) $args[] = "-p -B";
if($includeMultiOverlap) $args[] = "-O";
$tmp_dir = $parser->tempFolder();
$tmp_out = $parser->tempFile();
$args[] = "--tmpDir $tmp_dir";
$args[] = "-o $tmp_out";
$args[] = "$in";

//execute command
$parser->exec(get_path("feature_counts"), implode(" ", $args), true);

// copy output
$parser->exec("cp", "$tmp_out $out", true);
?>