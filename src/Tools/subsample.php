<?php

/**
	@page subsample
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("subsample", "Create a subsample with a specified number of reads.");

// mandatory arguments
$parser->addString("in_folder", "Input sample folder.", false);
$parser->addString("out_folder", "Output sample folder.", false);
$parser->addString("out_name", "Output sample name.", false);
$parser->addInt("reads", "Number of reads.", false);
$parser->addInt("seed", "Seed value for subsampling.", false, 4351);

// extract arguments
extract($parser->parse($argv));

// create output folder if necessary
if (! is_dir($out_folder))
{
	mkdir($out_folder);
}

// determine input files
$in_fastq = [];
$in_fastq["R1"] = glob("{$in_folder}/*_R1_001.fastq.gz");
$in_fastq["R2"] = glob("{$in_folder}/*_R2_001.fastq.gz");
// remove R1, R2 entries if no files found
foreach ($in_fastq as $key => $value)
{
	if (count($value) == 0)
	{
		unset($in_fastq[$key]);
	}
}

// run pipeline for R1, R2
foreach ($in_fastq as $key => $value)
{
	$out = "{$out_folder}/{$out_name}_{$key}_001.fastq.gz";

	$pipeline = [];

	$pipeline[] = ["zcat", implode(" ", $value)];
	$pipeline[] = ["seqtk", "sample -s {$seed} - {$reads}"];
	$pipeline[] = ["gzip -c -1", "> {$out}"];

	$parser->execPipeline($pipeline, "subsampling");
}