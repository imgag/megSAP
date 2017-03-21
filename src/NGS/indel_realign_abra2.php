<?php

/*
	@page indel_realign_abra2
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("indel_realign_abra2", "Perform indel realignment using ABRA2.");
$parser->addInfileArray("in",  "Input BAM file(s).", false);
$parser->addOutfileArray("out",  "Output BAM file(s).", false);
//optional
$parser->addString("build", "The genome build to use.", true, "hg19");
$parser->addInfile("roi", "Target region for realignment.", true, "");
$parser->addFlag("se", "Single-end input", true, "");
$parser->addInfile("gtf",  "GTF annotation file.", true, "");
$parser->addInfile("junctions",  "Junctions output file from STAR mapping.", true, "");
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
extract($parser->parse($argv));

//init
$local_data = get_path("local_data");
if (count($in)!=count($out))
{
	trigger_error("The number of input and output files must match!", E_USER_ERROR);
}

//indel realignment with ABRA
$tmp_folder = $parser->tempFolder("abra");
$tmp_unsorted = array();
for($i=0; $i<count($in); ++$i)
{
	$tmp_unsorted[] = $parser->tempFile(".bam");
}

$params = array();
$params[] = "--in ".implode(",", $in);
$params[] = "--out ".implode(",", $tmp_unsorted);
$params[] = "--threads $threads";
$params[] = "--tmpdir $tmp_folder";
$params[] = "--ref {$local_data}/{$build}.fa";
if ($gtf != "") $params[] = "--gtf $gtf";
if (isset($junctions)) $params[] = "--junctions $junctions";
if ($se) $params[] = "--single";
if (isset($roi)) $params[] = "--targets $roi";

$parser->exec(get_path("abra2"), implode(" ", $params), true);

//sort by position and create index file
$tmp_for_sorting = $parser->tempFile();
for($i=0; $i<count($out); ++$i)
{
	$parser->exec(get_path("samtools"), "sort -O bam -T $tmp_for_sorting -m 1G -@ ".min($threads, 4)." -o ".$out[$i]." ".$tmp_unsorted[$i], true);
	$parser->exec(get_path("ngs-bits")."BamIndex", "-in ".$out[$i], true);
}

?>