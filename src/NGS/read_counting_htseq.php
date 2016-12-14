<?php

/**
	@page read_counting_htseq
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("read_counting_htseq", "Perform read counting for aligned reads using htseq-count.");
$parser->addInfile("in",  "Input file in bam format.", false, "bam");
$parser->addString("prefix", "String to add to output files. Might include sub-directories.", false, NULL);

//optional parameters
$parser->addInt("threads", "Number of threads used for read counting", true, "4");
$parser->addString("method", "Read counting method", true, "intersection-nonempty");
$parser->addString("featureType", "Feature type used for mapping reads to features.", true, "exon");
$parser->addString("gtfAttribute", "GTF attribute used as feature ID.", true, "gene_id");
$parser->addString("gtfFile", "GTF File containing feature annotations used for read counting.", true, get_path("data_folder")."genomes/gtf/Homo_sapiens.GRCh37.75.renamed.gtf");
$parser->addFlag("stranded", "Specify whether a stranded protocol was used during library preparation. Default is non-stranded.");
$parser->addString("inputFormat", "Input file format. Options are 'sam' or 'bam'.", true, "bam"); //default in HTSeq is sam format, but it can be specified with "-f sam/bam"
$parser->addInt("minAQual", "Minimal alignment quality. Skip all reads with alignment quality lower than the given value.", true, 3);
$parser->addString("sortOrder", "Specify the sort order of the input sam/bam file. Input files have to be sorted either by name (name) or by alignment position (pos).", true, "pos");
$parser->addFlag("fixmate", "Run samtools fixmate before execution. bam file has to be sorted by names!");
$parser->addFlag("paired", "The data is paired-end. Only properly paired reads are used.");

extract($parser->parse($argv));

//extracting sub-directories and generating folder structure
$out=create_path($prefix);
$sampleName = $out[0];
$outdir = $out[1];

$parser->log("read_counting_htseq output directory=$outdir");

//build command
$arguments = array();

if($stranded) {
	$arguments[] = "-s yes";
} else {
	$arguments[] = "-s no";
}

$arguments[] = "-f $inputFormat";
$arguments[] = "-a $minAQual";
$arguments[] = "-t $featureType";
$arguments[] = "-i $gtfAttribute";
$arguments[] = "-m $method";

if($fixmate) {
	$parser->log("Starting fixmate operations on file: $in");
	if($sortOrder !== "name") {
		$parser->log("Sorting bam file by name");
		$in_sorted=$parser->tempFile(".nameSorted");
		$parser->exec(get_path("samtools"), "sort -@ $threads -n $in $in_sorted", true);
		$in_sorted .= ".bam";
		$parser->log("Starting fixmate on name sorted file: $in_sorted");
		$in_fixmate=$parser->tempFile(".nameSorted.fixmate.bam");
		$parser->exec(get_path("samtools"), "fixmate -r $in_sorted $in_fixmate", true);
		$parser->exec(get_path("samtools"), "index $in", true);
		$in="$in_fixmate";
	} else {
		$parser->log("Starting fixmate on name sorted file");
		$in_fixmate=$parser->tempFile(".fixmate.bam");
		$parser->exec(get_path("samtools"), "fixmate -r -p $in $in_fixmate", true);
		$parser->exec(get_path("samtools"), "index $in_fixmate", true);
		$in="$in_fixmate";
	}
	$arguments[] = "-r name";
} else {
	$arguments[] = "-r $sortOrder";
}

if($paired) {
	$parser->log("Extracting properly paired reads for counting");
	$in_pp=$parser->tempFile(".pp.bam");
	$parser->exec(get_path("samtools"), "view -b -f 2 $in > $in_pp", true);
	$in="$in_pp";
}

$arguments[] = "$in";
$arguments[] = "$gtfFile";
$arguments[] = "> ${outdir}${sampleName}.HTSeq.tsv"; //output file

//execute command
$parser->log("Starting read counting");
$parser->exec(get_path("htseq-count"), implode(" ", $arguments), true);
?>