<?php

/**
 * @page quant_salmon
 * 
 * @todo gene-level abundances --geneMap <gtf-file>
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("quant_salmon", "Transcript quantification with salmon.");

//mandatory arguments
$parser->addInfileArray("in1", "Input forward/R1 FASTQ file(s).", false);
$parser->addOutfile("out", "Abundance output file in tabular format.", false);

//optional arguments
$parser->addInfileArray("in2", "Input reverse/R2 FASTQ file(s).", true);
$parser->addEnum("library_type", "Library type specification, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addString("index", "Transcriptome index (salmon).", true, get_path("data_folder")."genomes/salmon/GRCh37_cDNA");
$parser->addInt("threads", "Number of threads.", true, 8);

//extract arguments
extract($parser->parse($argv));

//resolve library type parameter for salmon, see
//http://salmon.readthedocs.io/en/latest/library_type.html#fragment-library-types
$strandedness_paired = ["unstranded" => "IU",
	"reverse" => "ISR",
	"forward" => "ISF"
   ];
$strandedness_single = ["unstranded" => "U",
	"reverse" => "SR",
	"forward" => "SF"
   ];

//resolve fastq input parameters
if (isset($in2))
{
	$in = "--mates1 " . implode(" ", $in1) . " --mates2 " . implode(" ", $in2);
	$lib_type = $strandedness_paired[$library_type];
}
else
{
	$in = "--unmatedReads "  . implode(" ", $in1);
	$lib_type = $strandedness_single[$library_type];
}

$tmp_folder = $parser->tempFolder();

$args = ["--no-version-check",
	"quant",
	"--libType", $lib_type,
	$in,
	"--index", $index,
	"--threads", $threads,
	"--output", $tmp_folder
	];

$parser->exec(get_path("salmon"), implode(" ", $args), true);

$parser->log("salmon aux_info/meta_info.json:",
	file($tmp_folder . "/aux_info/meta_info.json"));

$parser->log("salmon lib_format_counts.json:",
	file($tmp_folder . "/lib_format_counts.json"));

$parser->copyFile($tmp_folder . "/quant.sf", $out);