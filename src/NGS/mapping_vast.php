<?php

/**
	@page mapping_vast
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("mapping_vast", "Align RNA-Seq data to VASTDB genome/junction libraries.");

// mandatory arguments
$parser->addString("sample_dir", "Sample folder with FASTQ files.", false);

// optional values
$parser->addString("out_dir", "Output directory, defaults to sample_dir/vast.", true, "");
$parser->addEnum("species", "Species/database to use.", true, ["Hsa", "Mmu"], "Mmu");
$parser->addString("db_dir", "Database directory.", true, "/mnt/users/ahmattj1/software/VASTDB");
$parser->addString("bowtie_prog", "path to bowtie program.", true, "/mnt/share/opt/bowtie-1.0.0/bowtie");
$parser->addInt("threads", "Number of threads for bowtie.", true, 4);

// extract arguments
extract($parser->parse($argv));

// resolve output directory
if ($out_dir === "")
{
	$out_dir = "{$sample_dir}/vast";
}
if (! is_dir($out_dir))
{
	mkdir($out_dir);
}

// extract sample name
$sample_name = preg_replace("/^Sample_/", "", basename($sample_dir));

// input FASTQ files
$in_fastq["R1"] = glob("{$sample_dir}/*_R1_001.fastq.gz");
$in_fastq["R2"] = glob("{$sample_dir}/*_R2_001.fastq.gz");

$merged_fastq = [];

// run concatenation pipeline for R1, R2
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
		$merged_fastq[] = $value[0];
	}
	// multiple FASTQ files, concatenation
	else
	{
		$out = $parser->tempFile("_{$key}.fastq.gz", "{$sample_name}_");

		$pipeline = [];

		$pipeline[] = ["zcat", implode(" ", $value)];
		$pipeline[] = ["gzip -c -1", "> {$out}"];

		$parser->execPipeline($pipeline, "FASTQ concatenation");

		$merged_fastq[] = $out;
	}
}

// run vast align
$vast_align_args = [
	"align",
	implode(" ", $merged_fastq),
	"--sp", $species,
	"--dbDir", $db_dir,
	"--bowtieProg", $bowtie_prog,
	"--expr",
	"--output", $out_dir,
	"--cores", $threads
];
$parser->exec("/mnt/SRV018/users/ahmattj1/software/vast-tools-1.3.0/vast-tools",
	implode(" ", $vast_align_args), true);

// run vast combine
$vast_combine_args = [
	"combine",
	"--output", $out_dir,
	"-sp", $species,
	"--dbDir", $db_dir,
	"-C"
];

// mouse: lift-over to mm10/GRCm38
if ($species === "Mmu")
{
	$vast_combine_args[] = "-a mm10";
}
$parser->exec("/mnt/SRV018/users/ahmattj1/software/vast-tools-1.3.0/vast-tools",
	implode(" ", $vast_combine_args), true);