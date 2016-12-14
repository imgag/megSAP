<?php

/*
	@page indel_realign_abra
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("indel_realign_abra", "Perform indel realignment using ABRA.");
$parser->addInfileArray("in",  "Input BAM file(s).", false);
$parser->addInfile("roi", "Target region for realignment", false);
$parser->addOutfileArray("out",  "Output BAM file(s).", false);
$parser->addFloat("mer",  "ABRA minimum edge pruning ratio parameter.", false);
//optional
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "hg19");
extract($parser->parse($argv));

//init
$local_data = get_path("local_data");
if (count($in)!=count($out))
{
	trigger_error("The number of input and output files must match!", E_USER_ERROR);
}

//create k-mer folder
$kmer_folder = get_path("data_folder")."/dbs/ABRA/";
if (!file_exists($kmer_folder) && !mkdir($kmer_folder, 0777, true))
{
	trigger_error("Could not create ABRA database folder '$kmer_folder'!", E_USER_ERROR);
}

//pre-calcalculate target region file with k-mer size (if not present)
$kmer_file = $kmer_folder.basename($roi, ".bed")."_".md5_file($roi).".bed";
if (!file_exists($kmer_file))
{
	$read_length = 100;
	$kmer_tmp = $parser->tempFile();
	$parser->exec(get_path("abra")." abra.KmerSizeEvaluator", "$read_length {$local_data}/{$build}.fa $kmer_tmp $threads $roi", true);
	copy2($kmer_tmp, $kmer_file);
}

//add BWA to path
putenv("PATH=".dirname(get_path("bwa")).":".getenv("PATH"));

//indel realignment with ABRA
$tmp_folder = $parser->tempFolder("abra");
$tmp_unsorted = array();
for($i=0; $i<count($in); ++$i)
{
	$tmp_unsorted[] = $parser->tempFile(".bam");
}
$parser->exec(get_path("abra")." abra.Abra", "--no-debug --in ".implode(",", $in)." --out ".implode(",", $tmp_unsorted)." --target-kmers $kmer_file --threads $threads --working $tmp_folder --ref {$local_data}/{$build}.fa --mer $mer -mad 5000", true);

//sort by position and create index file
$tmp_for_sorting = $parser->tempFile();
for($i=0; $i<count($out); ++$i)
{
	$parser->exec(get_path("samtools"), "sort -O bam -T $tmp_for_sorting -m 1G -@ ".min($threads, 4)." -o ".$out[$i]." ".$tmp_unsorted[$i], true);
	$parser->exec(get_path("ngs-bits")."BamIndex", "-in ".$out[$i], true);
}

?>