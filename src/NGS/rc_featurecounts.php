<?php

/**
	@page rc_featurecounts
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("rc_featurecounts", "Perform read counting for aligned reads using featureCount contained in the subread package.");
$parser->addInfile("in",  "BAM input file.", false, true);
$parser->addOutfile("out", "Raw count output TSV file.", false);

//optional parameters
$parser->addFlag("exon_level", "Perform exon-level read counting.");

$parser->addString("gtf_file", "GTF file containing feature annotations used for read counting.", true, get_path("data_folder")."/dbs/gene_annotations/GRCh38.gtf");
$parser->addString("feature_type", "Feature type used for mapping reads to features.", true, "exon");
$parser->addString("gtf_attribute", "Attribute used as feature ID.", true, "gene_id");
$parser->addEnum("library_type", "Specify the library type, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addFlag("single_end", "Single-end data.");
$parser->addInt("min_mapq", "Minimal mapping quality.", true, 3);
$parser->addFlag("overlap", "Count reads multiple times if they overlap more than one feature.");
$parser->addString("qc_file", "If specified, write expression/assignment QC file.", true, "");
$parser->addFlag("qc_skip_strandedness", "Skip strandedness/library compatibility QC calculations.");
$parser->addFlag("ignore_dup", "Do not count alignments flagged as duplicates.");
$parser->addOutfile("out_report", "Output report file in BAM format.", true);
$parser->addInt("threads", "Number of threads used for read counting", true, 4);
extract($parser->parse($argv));

$strandedness = [
	"unstranded" => 0,
	"reverse" => 2,
	"forward" => 1
];

$reverse_library_type = [
	"unstranded" => "unstranded",
	"reverse" => "forward",
	"forward" => "reverse"
];

$tmp_dir = $parser->tempFolder();
$tmp_dir2 = $parser->tempFolder();
$tmp_out = "{$tmp_dir2}/counts";

//build arguments array
$args = [
	"-a", $gtf_file,
	"-t", $feature_type,
	"-g", $gtf_attribute,
	"-Q", $min_mapq,
	"-T", $threads,
	"-s", $strandedness[$library_type],
	"--tmpDir", $tmp_dir,
	"-o", $tmp_out
];

if (!$single_end)
{
	$args[] = "-p -B --countReadPairs";
}
if ($overlap || $exon_level)
{
	$args[] = "-O";
}
if ($exon_level)
{
	$args[] = "-f";
}
if ($ignore_dup)
{
	$args[] = "--ignoreDup";
}
if (isset($out_report))
{
	$args[] = "-R BAM";
}
$args[] = $in;

// execute command
$parser->exec(get_path("feature_counts"), implode(" ", $args), true);

//write summary into the tool log
$parser->log("featureCounts summary", file("{$tmp_out}.summary"));

// read counts TSV
$parser->copyFile($tmp_out, $out);
$parser->copyFile("{$tmp_out}.summary", "${out}.summary");
// report BAM file
if (isset($out_report))
{
	$in_basename = basename($in);
	$report_bam = "{$tmp_dir2}/{$in_basename}.featureCounts.bam";
	$sort_tmp = $parser->tempFile(".bam");
	$parser->exec(get_path("samtools"), "sort -T {$sort_tmp} -m 1G -@ 4 -o {$out_report} {$report_bam}", true);
	$parser->indexBam($out_report, $threads);
}

if ($qc_file !== "")
{
	$qc = new Matrix();


	// QC for all libraries, based on counts.summary
	$summary = array_column(load_tsv("${tmp_out}.summary"), 1, 0);

	// number of usable fragments
	$usable_fragments = $summary["Assigned"] +
		$summary["Unassigned_Ambiguity"] +
		$summary["Unassigned_NoFeatures"];
	$qc->addRow([
		"usable fragments number",
		$usable_fragments
	]);

	// % assigned reads
	$uniq_ass_percentage = 100.0 * $summary["Assigned"] / $usable_fragments;
	$qc->addRow([
		"unique assignment percentage",
		$uniq_ass_percentage
	]);

	// % ambiguous reads
	$qc->addRow([
		"ambiguous assignment percentage",
		100.0 * $summary["Unassigned_Ambiguity"] / $usable_fragments
	]);

	// % no feature
	$qc->addRow([
		"no assignment percentage",
		100.0 * $summary["Unassigned_NoFeatures"] / $usable_fragments
	]);

	// library compatibility / strand bias

	// for unstranded libraries
	if ($library_type === "unstranded" && !$qc_skip_strandedness)
	{
		$uniqe_assignment_rates = [];
		foreach (["forward", "reverse"] as $lib_type)
		{
			$qc_tmp_out = $parser->tempFile();
			$qc_tmp_qc = $parser->tempFile();

			$qc_args = [
				"-in", $in,
				"-out", $qc_tmp_out,
				"-qc_file", $qc_tmp_qc,
				"-library_type", $lib_type,
				"-gtf_file", $gtf_file,
				"-gtf_attribute", $gtf_attribute,
				"-feature_type", $feature_type,
				"-threads", $threads,
				"-qc_skip_strandedness"
			];
			if ($single_end)
			{
				$qc_args[] = "-single_end";
			}
			if ($overlap)
			{
				$qc_args[] = "-overlap";
			}
			if ($exon_level)
			{
				$qc_args[] = "-exon_level";
			}
			if ($ignore_dup)
			{
				$qc_args[] = "-ignore_dup";
			}

			$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $qc_args));


			$uniqe_assignment_rates[] = array_column(load_tsv($qc_tmp_qc), 1, 0)["unique assignment percentage"];
		}

		$qc->addRow(["forward library compatibility", $uniqe_assignment_rates[0]]);
		$qc->addRow(["reverse library compatibility", $uniqe_assignment_rates[1]]);
	}
	// for stranded library
	// count with reversed library
	else if (!$qc_skip_strandedness)
	{
		$rev_lib_type = $reverse_library_type[$library_type];

		$qc_tmp_out = $parser->tempFile();
		$qc_tmp_qc = $parser->tempFile();

		$qc_args = [
			"-in", $in,
			"-out", $qc_tmp_out,
			"-qc_file", $qc_tmp_qc,
			"-library_type", $rev_lib_type,
			"-gtf_file", $gtf_file,
			"-gtf_attribute", $gtf_attribute,
			"-feature_type", $feature_type,
			"-threads", $threads,
			"-qc_skip_strandedness"
		];
		if ($single_end)
		{
			$qc_args[] = "-single_end";
		}
		if ($overlap)
		{
			$qc_args[] = "-overlap";
		}
		if ($ignore_dup)
		{
			$qc_args[] = "-ignore_dup";
		}

		$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $qc_args));


		$qc->addRow([$library_type . " library compatibility",
			$uniq_ass_percentage]);
		$qc->addRow([$rev_lib_type . " library compatibility",
			array_column(load_tsv($qc_tmp_qc), 1, 0)["unique assignment percentage"]]);
	}
	
	// number of expressed genes
	$cutoff_values = [1, 3, 10];
	foreach ($cutoff_values as $cutoff)
	{
		$pipeline = [];
		$pipeline[] = ["tail", "-n+3 {$tmp_out}"];
		$pipeline[] = [get_path("ngs-bits") . "TsvFilter", "-numeric -filter '7 >= {$cutoff}'"];
		$pipeline[] = ["wc", "-l"];
		$ret = $parser->execPipeline($pipeline, "count-genes");
		
		$num_genes = intval($ret[0][0]) - 1;
		$qc->addRow(["number of genes with at least {$cutoff} read(s)", $num_genes]);
	}

	$qc->setHeaders(["key", "value"]);
	$qc->toTSV($qc_file);
}