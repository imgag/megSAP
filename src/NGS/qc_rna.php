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
$parser->addInfile("system",  "Processing system INI file (determined from file name by default).", true);

//extract arguments
extract($parser->parse($argv));

$processed_sample = basename($in, ".bam");
// extract processing system information
$sys = load_system($system, $processed_sample);
if (!isset($sys['build']) || $sys['build']=="") {
	trigger_error("Build not specified!", E_USER_ERROR);
} else {
	$build = $sys['build'];
}

$gtf = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

$bam = $in;

if (!is_dir($out))
{
	mkdir($out);
}

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
//$results->setHeaders(array("key", "value"));
$results->addRow(array("input", $reads_input));
$results->addRow(array("unique", $reads_uniq));
$results->addRow(array("multi", $reads_multimapped));
$results->addRow(array("unmapped", $reads_unmapped));

$results->toTSV("{$out}/mapping.tsv");

$parser->log(">>> Calculating chromosome distribution using samtools index");
$parser->exec(get_path("samtools"), "idxstats {$bam} > {$out}/idxstats.tsv", true);

$parser->log(">>> Calculating gene expression and biotypes using rc_ scripts");
$parser->execTool("NGS/rc_featurecounts.php", "-in {$bam} -out {$out}/featurecounts_raw.tsv -keep_summary -library_type reverse -gtf_file {$gtf}");
$parser->execTool("NGS/rc_normalize.php", "-in {$out}/featurecounts_raw.tsv -out {$out}/gene_expression.tsv -method raw");
$parser->execTool("NGS/rc_annotate.php", "-in {$out}/gene_expression.tsv -out {$out}/gene_expression.tsv -annotationId gene_biotype -gtfFile {$gtf}");
$parser->exec("rm", "{$out}/featurecounts_raw.tsv", true);
$parser->exec("mv", "{$out}/featurecounts_raw_summary.tsv {$out}/featurecounts_summary.tsv", true);