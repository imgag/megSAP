<?php 

/** 
	@page merge_fastq
*/


require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("merge_fastq", "Merges FASTQ files if a sample was distributed over several runs/lanes/MIDs or processings. Merges on a per sample basis.");
$parser->addInfile("in",  "Input sample folder (e.g. 'XYZ/Sample_GS99999_99'). Only fastq files of the sample within the folder name will be merged.", false);
// optional
$parser->addString("out",  "Output folder. Default is same as in. If specified otherwise data will only be copied (and not moved).", true, "default");
extract($parser->parse($argv));

//check/create input, output and backup folder
$in = realpath($in)."/";
if(!is_dir($in)) trigger_error("Cannot access $in", E_USER_ERROR);
if($out=="default")	$out = $in;
$out = realpath($out)."/";
if(!is_dir($out)) trigger_error("Cannot access $out", E_USER_ERROR);

$backup = "+original/";
if (!file_exists($out.$backup)) mkdir($out.$backup);

//determine sample name from FASTQ files
$s_id = "";
$fastq_files = array_merge(glob($in."*.fastq.gz"), glob($in.$backup."*.fastq.gz"));
foreach($fastq_files as $file)
{
	$file = basename($file);
	list($tmp) = explode("_", $file);
	if ($s_id=="") $s_id = $tmp;
	if ($s_id!=$tmp) trigger_error("FASTQ files from different samples ($s_id and $tmp) cannot be merged!", E_USER_ERROR);
}

//identify fastq files and move / copy original files to backup folder
$files = array();
$tmp_files = array_unique(array_merge(glob($in.$s_id."*.fastq.gz"), glob($in.$backup.$s_id."*.fastq.gz"), glob($out.$backup.$s_id."*.fastq.gz")));	//array_unique required since $in and $out can be identical
foreach($tmp_files as $tmp_index => $tmp_file)
{
	if(contains($tmp_file, "merged") || contains($tmp_file, "trimmed"))	continue;
	if($in==$out && !contains($tmp_file, $backup))
	{
		$parser->moveFile($tmp_file, $in.$backup);
	}
	else if($in!=$out)
	{
		$parser->copyFile($tmp_file, $out.$backup);
	}
	$files[] = basename($tmp_file);
}

//sort files into groups according to read number (1=forward, 2=reverse)
$groups = array();
foreach ($files as $file)
{
	$parts = explode("_", $file);
	$name = $parts[0];
	$read = $parts[count($parts)-2];
	$group = $name."_merged_".$read."_001.fastq.gz";
	if (!isset($groups[$group])) $groups[$group] = array();
	$groups[$group][] = $file;
}

//merge files
foreach($groups as $name => $files)
{
	$param = "";
	foreach($files as $file)
	{
		$param .= $out.$backup.$file." ";
	}

	$parser->exec("zcat", $param." | gzip -1 > ".$out."/".$name, true);
}