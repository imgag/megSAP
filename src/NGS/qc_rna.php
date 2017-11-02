<?php

/**
	@page qc_rna
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("qc_rna", "Calculate several RNA-seq QC parameters.");

//mandatory arguments
$parser->addInfile("in", "input BAM file", false, true);
$parser->addString("out", "output folder", false);

//optional arguments
$parser->addInfile("system",  "Processing system INI file (determined from BAM file name by default).", true);
$parser->addString("steps", "QC steps to perform\nbam: alignment metrics; expr: expression metrics; dup: duplication rate analysis", true, "bam,expr,dup");
$parser->addFlag("se", "flag indicating single-end sample");
$parser->addEnum("library_type", "Specify the library type, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);

//extract arguments
extract($parser->parse($argv));
$steps = explode(",", $steps);

// extract processing system information
$processed_sample = basename($in, ".bam");
$sys = load_system($system, $processed_sample);
if (!isset($sys['build']) || $sys['build'] === "")
{
	trigger_error("Build not specified!", E_USER_ERROR);
}
else
{
	$build = $sys['build'];
}

$bam = $in;
$gtf = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";
$paired = $se ? "no" : "yes";
$strandedness = [
	"unstranded" => "0",
	"reverse" => "2",
	"forward" => "1"
   ];
$lib_type = $strandedness[$library_type];

if (!is_dir($out))
{
	mkdir($out);
}

// BAM/alignment metrics
if (in_array("bam", $steps))
{
	$parser->log(">>> Calculating BAM metrics using samtools index/view...");

	//samtools: unique alignments (= STAR uniquely mapped reads)
	$reads_uniq = $parser->exec(get_path("samtools"), "view -q 255 -c {$bam}", true);
	$reads_uniq = intval($reads_uniq[0][0]);

	//samtools: secondary alignments
	$reads_secondary = $parser->exec(get_path("samtools"), "view -f0x100 -c {$bam}", true);
	$reads_secondary = intval($reads_secondary[0][0]);

	//samtools: all alignments
	$reads_all = $parser->exec(get_path("samtools"), "idxstats {$bam}| cut -f3 | paste -s -d'+' | bc", true);
	$reads_all = intval($reads_all[0][0]);

	//samtools: unmapped reads
	$reads_unmapped = $parser->exec(get_path("samtools"), "idxstats {$bam}| cut -f4 | paste -s -d'+' | bc", true);
	$reads_unmapped = intval($reads_unmapped[0][0]);

	$reads_multimapped = $reads_all - $reads_secondary - $reads_uniq;
	$reads_input = $reads_all - $reads_secondary + $reads_unmapped;

	$results = new Matrix();
	$results->addRow(array("input", $reads_input));
	$results->addRow(array("unique", $reads_uniq));
	$results->addRow(array("multi", $reads_multimapped));
	$results->addRow(array("unmapped", $reads_unmapped));

	$results->toTSV("{$out}/mapping.tsv");

	$parser->log(">>> Calculating chromosome distribution using samtools index");
	$parser->exec(get_path("samtools"), "idxstats {$bam} > {$out}/idxstats.tsv", true);
}

// expression/read count metrics
if (in_array("expr", $steps))
{
	$parser->log(">>> Calculating gene expression and biotypes using rc_ scripts");
	
	if ($library_type === "unstranded")
	{
		// sense+antisense expression
		$parser->execTool("NGS/rc_featurecounts.php", "-in {$bam} -out {$out}/featurecounts_raw.tsv -keep_summary -library_type unstranded -gtf_file {$gtf}");
	}
	else
	{		
		// sense expression
		$parser->execTool("NGS/rc_featurecounts.php", "-in {$bam} -out {$out}/featurecounts_raw.tsv -keep_summary -library_type reverse -gtf_file {$gtf}");

		// antisense expression
		$parser->execTool("NGS/rc_featurecounts.php", "-in {$bam} -out {$out}/featurecounts_raw_antisense.tsv -keep_summary -library_type forward -gtf_file {$gtf}");
	}
	
	
	// create tsv from sense expression
	$parser->execTool("NGS/rc_normalize.php", "-in {$out}/featurecounts_raw.tsv -out {$out}/gene_expression.tsv -method raw");
	
	// annotate sense expression with biotype
	$parser->execTool("NGS/rc_annotate.php", "-in {$out}/gene_expression.tsv -out {$out}/gene_expression.tsv -annotationId gene_biotype -gtfFile {$gtf}");
	
	// remove raw sense expression
	$parser->exec("rm", "{$out}/featurecounts_raw.tsv", true);
	
	// keep summary of sense expression
	$parser->exec("mv", "{$out}/featurecounts_raw_summary.tsv {$out}/featurecounts_summary.tsv", true);
}

// duplication rate analysis
if (in_array("dup", $steps))
{
	$script = <<<'EOT'
options(bitmapType='cairo')
library(dupRadar)
args <- commandArgs(TRUE)

bam <- args[1]

gtf <- gsub("gtf=","",args[2])

## no=0|yes=1|reverse=2
stranded <- as.integer(gsub("stranded=","",args[3]))

## is a paired end experiment (yes/no)
paired   <- gsub("paired=","",args[4])

## output directory
outdir   <- gsub("outdir=","",args[5])

## number of threads to be used
threads  <- as.integer(gsub("threads=","",args[6]))


dm <- analyzeDuprates(bam, gtf, stranded, (paired == "yes"), threads)

## duprate vs. expression smooth scatter
png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_duprate.png"),
    width=1000, height=1000)
duprateExpDensPlot(dm, main = basename(bam))
silent <- dev.off()

## fit intercept and slope numeric values
fit <- duprateExpFit(DupMat=dm)
cat(sprintf("fit.intercept\t%f\nfit.slope\t%f\n", fit$intercept, fit$slope))
EOT;
	
	$scriptfile = $parser->tempFile("_dupradar.R");
	file_put_contents($scriptfile, $script);
	$script_args = [
		$scriptfile,
		$bam,
		$gtf,
		"stranded={$lib_type}",
		"paired={$paired}",
		"outdir={$out}",
		"threads={$threads}",
		" > {$out}/dupradar.tsv"
	];
	$parser->exec(get_path("rscript"), implode(" ", $script_args), true);
}