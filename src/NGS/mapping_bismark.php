<?php

/**
	@page mapping_bismark
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("mapping_bismark", "Bisulfit read mapping using bismark.");

// mandatory arguments
$parser->addString("folder", "Sample folder containing FASTQ files.", false);

// optional arguments
$parser->addString("out_folder", "Folder where results are stored, defaults to sample folder.", true, "default");
$parser->addInfile("system",  "Processing system INI file (determined from folder name by default).", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);

extract($parser->parse($argv));

// resolve output folder
if ($out_folder=="default")
{
	$out_folder = $folder;
}

// input FASTQ files
$in_for = glob($out_folder."/*_R1_001.fastq.gz");
$in_rev = glob($out_folder."/*_R2_001.fastq.gz");

if ((count($in_for) == 0) && (count($in_rev) == 0)) {
	trigger_error("No FASTQ files found!", E_USER_ERROR);
}
if ((count($in_for) != count($in_rev)) && count($in_rev)>0) {
	trigger_error("Mismatching number of R1 and R2 files!", E_USER_ERROR);
}

$paired = count($in_rev) > 0;

// output file names
// use sample id, if possible
$p = basename(realpath($out_folder));
$name = "";
if (preg_match('/^Sample_(.+)/', $p, $matches)) {
	$name = $matches[1];
}

$sys = load_system($system, $name);

$in_for_s = implode(" ", $in_for);
$in_rev_s = implode(" ", $in_rev);

// perform adapter trimming, read QC and merging
if ($sys["adapter1_p5"]=="" && $sys["adapter2_p7"]=="")
{
	trigger_error("No forward and/or reverse adapter sequence given!\nForward: ".$sys["adapter1_p5"]."\nReverse: ".$sys["adapter2_p7"], E_USER_ERROR);
}
$trimmed1 = $parser->tempFile("_trimmed1.fastq.gz");
$trimmed2 = $parser->tempFile("_trimmed2.fastq.gz");
$basename = $out_folder."/".$name;
$stafile1 = $basename."_stats_fastq.qcML";
$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 $in_for_s -in2 $in_rev_s -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -threads $threads", true);

// run bismark
$bismark_tmp = $parser->tempFolder();
$tmp_out_folder = $parser->tempFolder();

$bismark_params = array();
$bismark_params[] = get_path("data_folder")."/genomes/bismark/".$sys['build'];
$bismark_params[] = "-1 $trimmed1 -2 $trimmed2";
$bismark_params[] = "-N 1 --unmapped";
$bismark_params[] = "--bowtie2 --path_to_bowtie ".dirname(get_path("bowtie2"));
$bismark_params[] = "--rg_tag TRUE --rg_id SAMPLE --rg_sample $name";
$bismark_params[] = "--basename $name";
$bismark_params[] = "-p $threads";
$bismark_params[] = "--temp_dir $bismark_tmp";
$bismark_params[] = "-o $tmp_out_folder";

$parser->exec(get_path("bismark"), implode(" ", $bismark_params), true);

// bismark output files
$bismark_bam = $tmp_out_folder."/".$name."_pe.bam";
$parser->exec(dirname(get_path("bismark"))."/bismark_methylation_extractor", "-p --no_overlap --multicore $threads --gzip --output $out_folder --bedGraph --cytosine_report $bismark_bam", true);
$parser->exec(dirname(get_path("bismark"))."/deduplicate_bismark", "-p --bam $bismark_bam", true);
$parser->exec(dirname(get_path("bismark"))."/bismark2report", "--dir $tmp_out_folder --alignment_report {$tmp_out_folder}/{$name}_PE_report.txt", true);

$parser->moveFile($tmp_out_folder."/".$name."_PE_report.html", $out_folder."/".$name."_bismark_report.html");
$parser->moveFile($tmp_out_folder."/".$name."_pe.bismark.cov.gz", $out_folder."/".$name."_CpG.tsv.gz");

//sort and index BAM
$out = $out_folder."/".$name.".bam";
$parser->sortBam($bismark_bam, $out, $threads);
$parser->indexBam($out, $threads);

//TODO MappingQC

?>
