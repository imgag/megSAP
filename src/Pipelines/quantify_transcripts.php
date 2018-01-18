<?php

/**
	@page quantify_transcripts
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("quantify_transcripts", "RNA transcript quantification pipeline using salmon.");

//mandatory
$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);

//optional
$parser->addInfile("system", "Processing system INI file (determined from NGSD via the 'name' by default).", true);
$parser->addEnum("library_type", "Specify the library type, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
$parser->addString("out", "Output file, defaults to <name>_transcript_quant.tsv.", true, "");
$parser->addString("out_raw", "Raw output file, defaults to <name>_transcript_quant.sf.", true, "");

extract($parser->parse($argv));

//resolve out_folder
if ($out === "")
{
	$out = "{$folder}/{$name}_transcript_quant.tsv";
}

if ($out_raw === "")
{
	$out_raw = "{$folder}/{$name}_transcript_quant.sf";
}

//log server name
list($server) = exec2("hostname -f");
$user = exec('whoami');
$parser->log("Executed on server: ".implode(" ", $server)." as ".$user);

//init
$sys = load_system($system, $name);
$build = $sys['build'];

//determine transcriptome index from build
$index = get_path("data_folder")."/genomes/salmon/{$build}_cDNA/";

//determine gtf from build
$gtfFile = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

//find FASTQ files
$in_for = glob($folder."/*_R1_001.fastq.gz");
$in_rev = glob($folder."/*_R2_001.fastq.gz");

if ((count($in_for) == 0) && (count($in_rev) == 0))
{
	trigger_error("No FASTQ files found!", E_USER_ERROR);
}

$paired = (count($in_rev) != 0);

//check FASTQ quality encoding
$files = array_merge($in_for, $in_rev);
foreach ($files as $file)
{
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."FastqFormat", "-in {$file}", true);
	if (!contains($stdout[2], "Sanger"))
	{
		trigger_error("Input file '{$file}' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
	}
}

//check that adapters are specified
if ($sys["adapter1_p5"]=="" || $sys["adapter2_p7"]=="")
{
	trigger_error("No forward and/or reverse adapter sequence given!\nForward: ".$sys["adapter1_p5"]."\nReverse: ".$sys["adapter2_p7"], E_USER_ERROR);
}

//adapter trimming + QC (SeqPurge for paired-end, ReadQC+skewer for single-end)
$qc_fastq = "{$folder}/{$name}_stats_fastq.qcML";
if ($paired)
{
	$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
	$fastq_trimmed2 = $parser->tempFile("_trimmed.fastq.gz");
	$seqpurge_params = array(
		"-in1", implode(" ", $in_for),
		"-in2", implode(" ", $in_rev),
		"-out1 {$fastq_trimmed1}",
		"-out2 {$fastq_trimmed2}",
		"-a1", $sys["adapter1_p5"],
		"-a2", $sys["adapter2_p7"],
		"-qc", $qc_fastq,
		"-threads", min($threads, 2), // use at most 2 threads
		"-qcut 0"
		);
	$parser->exec(get_path("ngs-bits")."SeqPurge", implode(" ", $seqpurge_params), true);
}
else
{
	$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 ".implode(" ", $in_for)." -out $qc_fastq", true);

	$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");

	$skewer_params = ["-x", $sys["adapter1_p5"],
		"-y", $sys["adapter2_p7"],
		"-m any",
		"--threads", min($threads, 2), // use at most 2 threads
		"--stdout",
		"-"
		];

	$pipline = [];
	$pipeline[] = array("zcat", implode(" ", $in_for));
	$pipeline[] = array(get_path("skewer"), implode(" ", $skewer_params));
	$pipeline[] = array("gzip", "-1 > {$fastq_trimmed1}");
	$parser->execPipeline($pipeline, "skewer");
}

//quantification
$args_quant = ["-threads", $threads,
	"-in1", $fastq_trimmed1,
	"-index", $index,
	"-library_type", $library_type,
	"-out", $out_raw
   ];
if ($paired)
{
	$args_quant[] = "-in2 $fastq_trimmed2";
}

$parser->execTool("NGS/quant_salmon.php", implode(" ", $args_quant));

$tmp_quant = $parser->tempFile("_quant.tsv");
$quant = Matrix::fromTSV($out_raw);
$quant->removeRow(0);
$quant->removeCol(1);
$quant->removeCol(2);
$quant->setHeaders(["transcript_id", "tpm", "estimated_count"]);
$quant->toTSV($tmp_quant);

//annotation
$args_annotate = ["-in", $tmp_quant,
	"-out", $out,
	"-gtfFile", $gtfFile,
	"-column_name", "transcript_id",
	"-keyId", "transcript_id",
	"-annotationIds", "transcript_name,gene_biotype",
	"-ignore_version_suffix"
	];
$parser->execTool("NGS/rc_annotate.php", implode(" ", $args_annotate));