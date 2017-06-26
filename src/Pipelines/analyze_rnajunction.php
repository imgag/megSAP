<?php

/**
	@page analyze_rnajunction
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("analyze_rnajunction", "RNA junction detection pipeline using subjunc.");

$parser->addInfileArray("in_for", "Forward reads in fastq.gz file(s).", false);
$parser->addString("out_folder", "Output folder.", false);
$parser->addString("out_name", "Output base file name, typically the processed sample ID (e.g. 'GS120001_01').", false, NULL);

//optional
$parser->addInfileArray("in_rev", "Reverse reads fastq.gz file(s) for paired end alignment.", true);
$parser->addInfile("system", "Processing system INI file (determined from 'out_name' by default).", true);
$parser->addString("genome", "Subread genome index, by default determined from system/build.", true, "");
$parser->addInfile("gtfFile", "GTF file containing feature annotations (for read counting).", true, "");

$parser->addInt("threads", "The maximum number of threads used.", true, 4);

extract($parser->parse($argv));

//log server name
list($server) = exec2("hostname -f");
$user = exec('whoami');
$parser->log("Executed on server: ".implode(" ", $server)." as ".$user);

//output folder
if (!is_dir($out_folder) && !mkdir($out_folder, 0777, true))
{
	trigger_error("Could not create output folder!", E_USER_ERROR);
}

//init
$prefix = $out_folder."/".$out_name;
$sys = load_system($system, $out_name);
$build = $sys['build'];
$target_file = $sys['target_file'];
$paired = isset($in_rev);

//determine gtf from build
if (!isset($gtfFile))
{
	$gtfFile = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";
}
//determine genome from build
if ($genome == "")
{
	$genome = get_path("data_folder")."/genomes/subread/{$build}";
}

//mapping and QC
$final_bam = $prefix.".bam";
$qc_fastq = $prefix."_stats_fastq.qcML";
$qc_map = $prefix."_stats_map.qcML";

//check FASTQ quality encoding
$files = $paired ? array_merge($in_for, $in_rev) : $in_for;
foreach($files as $file)
{
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."FastqFormat", "-in $file", true);
	if (!contains($stdout[2], "Sanger"))
	{
		trigger_error("Input file '$file' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
	}
}

//check that adapters are specified
if ($sys["adapter1_p5"]=="" || $sys["adapter2_p7"]=="")
{
	trigger_error("No forward and/or reverse adapter sequence given!\nForward: ".$sys["adapter1_p5"]."\nReverse: ".$sys["adapter2_p7"], E_USER_ERROR);
}

//adapter trimming + QC (SeqPurge for paired-end, Skewer/ReadQC for single-end)
if($paired)
{
	$fastq_trimmed1 = $parser->tempFile("_trimmed_R1.fastq.gz");
	$fastq_trimmed2 = $parser->tempFile("_trimmed_R2.fastq.gz");
	$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $fastq_trimmed1 -out2 $fastq_trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $qc_fastq -threads ".$threads." -qcut 0", true);
}
else
{
	$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 ".implode(" ", $in_for)." -out $qc_fastq", true);

	$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
	$skewer_stderr = $parser->tempFile("_skewer_stderr");
	$parser->exec("zcat", implode(" ", $in_for)." | ".get_path("skewer")." -x ".$sys["adapter1_p5"]." -y ".$sys["adapter2_p7"]." -m any --threads $threads --quiet --stdout -"." 2> $skewer_stderr | gzip -1 > $fastq_trimmed1", true);
	$parser->log("skewer log", file($skewer_stderr));
}

//mapping
$args = array("-out {$final_bam}",
	"-in1 {$fastq_trimmed1}",
	"-gtf {$gtfFile}",
	"-genome {$genome}",
	"-threads {$threads}");
if($paired) $args[] = "-in2 {$fastq_trimmed2}";

$parser->execTool("NGS/mapping_subjunc.php", implode(" ", $args));

//mapping QC
if (isset($target_file) && $target_file != "") {
	$mappingqc_target = "-roi {$target_file}";
} else {
	$mappingqc_target = "-rna";
}
$parser->exec(get_path("ngs-bits")."MappingQC", "-in {$final_bam} -out {$qc_map} {$mappingqc_target}", true);

?>