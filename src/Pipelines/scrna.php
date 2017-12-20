<?php

/**
 * @page scrna
 *
 * @TODO consider mapping and counting all R2 and apply cell-filtering at end
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("scrna", "Single-Cell RNA analysis pipeline for SureCell 3' WTA.");

// mandatory arguments
$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name.", false);

// optional arguments
$parser->addString("out_folder", "Output folder, defaults to <folder>/scrna.", true, "");
$parser->addInfile("system", "Processing system file (determined from NGSD via 'name' by default).", true);
$parser->addInt("cell_number", "Set defined number of cells.", true, -1);
$parser->addInt("expected_cell_number", "Set expected number of cells.", true, -1);
$steps_all = [ "whitelist", "extract", "mapping" ];
$parser->addString("steps", "Comma-separated list of steps to perform:\n" .
	"whitelist=cell barcode determination,\n" .
	"extract=barcode extraction and FASTQ generation,\n" .
	"mapping=mapping and cell-specific expression analysis",
	true, implode(",", $steps_all));
$parser->addInt("threads", "Number of parallel threads for STAR.", true, 4);

// extract arguments
extract($parser->parse($argv));

$steps = explode(",", $steps);

// resolve output directory
if ($out_folder === "")
{
	$out_folder = "{$folder}/scrna";
}
if (! is_dir($out_folder))
{
	mkdir($out_folder);
}

$bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNN";

// processing system
$sys = load_system($system, $name);
$build = $sys['build'];
$target_file = $sys['target_file'];

// determine genome from build
$genome = get_path("data_folder")."/genomes/STAR/{$build}/";

// determine gtf from build
$gtfFile = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

if (in_array("whitelist", $steps) || in_array("extract", $steps))
{
	// find FASTQ files
	$in_fastq["R1"] = glob("{$folder}/*_R1_001.fastq.gz");
	
	if (in_array("extract", $steps))
	{
		$in_fastq["R2"] = glob("{$folder}/*_R2_001.fastq.gz");
	}

	// merge R1, R2 files if necessary
	// run concatenation pipeline for R1, R2
	$in_tmp = [];
	foreach ($in_fastq as $key => $value)
	{
		// no FASTQ files
		if (count($value) === 0)
		{
			continue;
		}
		// one FASTQ file, use it as-is
		else if (count($value) === 1)
		{
			$in_tmp[$key] = $value[0];
		}
		// multiple FASTQ files, concatenation
		else
		{
			$out = $parser->tempFile("_{$key}.fastq.gz", "{$name}_");

			$pipeline = [];

			$pipeline[] = ["zcat", implode(" ", $value)];
			$pipeline[] = ["gzip -c -1", "> {$out}"];

			$parser->execPipeline($pipeline, "FASTQ concatenation");

			$in_tmp[$key] = $out;
		}
	}
	$in_R1_tmp = $in_tmp["R1"];
	$in_R2_tmp = $in_tmp["R2"];
}

if (in_array("whitelist", $steps))
{
	// (1) UMI-counts whitelist
	// pipe R1 (cell barcode + UMI) to umi_tools whitelist
	$out_whitelist =  "{$out_folder}/whitelist.txt";
	$out_plots = "{$out_folder}/expect_whitelist";

	$args_whitelist = [
		"whitelist",
		"--bc-pattern", $bc_pattern,
		"--plot-prefix", $out_plots,
		"--stdin", $in_R1_tmp,
		"--stdout", $out_whitelist,
		"--log2stderr",
	];
	if ($expected_cell_number > 0)
	{
		$args_whitelist[] = "--expect-cells";
		$args_whitelist[] = $expected_cell_number;
	}
	if ($cell_number > 0)
	{
		$args_whitelist[] = "--set-cell-number";
		$args_whitelist[] = $cell_number;
	}
	$parser->exec(get_path("umi_tools"), implode(" ", $args_whitelist), true);
}

$out_extracted = "{$out_folder}/{$name}_extracted_R2_001.fastq.gz";
if (in_array("extract", $steps))
{
	// (2) UMI-counts extract
	// extract barcodes and append to read names
	$args_extract = [
		"extract",
		"--bc-pattern", $bc_pattern,
		"--stdin", $in_R1_tmp,
		"--stdout", $out_extracted,
		"--read2-in", $in_R2_tmp,
		"--read2-stdout"
	];
	$parser->exec(get_path("umi_tools"), implode(" ", $args_extract), true);
}

if (in_array("mapping", $steps))
{
	// (3) RNA-seq mapping
	// using $out_extracted as input for mapping_star
	$bam_tmp = $parser->tempFile(".bam", $name);
	$args_star = [
		"-in1", $out_extracted,
		"-out", $bam_tmp,
		"-genome", $genome,
		"-skip_dedup",
		"-threads", $threads
	];
	$parser->execTool("NGS/mapping_star.php", implode(" ", $args_star));

	// (4) featureCounts with BAM output
	$bam_assigned = "{$out_folder}/{$name}_assigned.bam";
	$counts_sum ="{$out_folder}/{$name}_counts_sum.tsv";
	$qc_counts_sum ="{$out_folder}/{$name}_stats_counts_sum.tsv";

	$args_featurecounts = [
		"-in", $bam_tmp,
		"-out", $counts_sum,
		"-out_report", $bam_assigned,
		"-gtf_file", $gtfFile,
		"-library_type", "unstranded",
		"-single_end",
		"-qc_file", $qc_counts_sum
	];
	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args_featurecounts));

	// (5) count reads per gene and cell
	$out_counts = "{$out_folder}/{$name}_counts_cell.tsv";
	$args_count = [
		"count",
		"--per-gene",
		"--gene-tag", "XT",
		"--per-cell",
		"-I", $bam_assigned,
		"-S", $out_counts
	];
	$parser->exec(get_path("umi_tools"), implode(" ", $args_count), true);
}