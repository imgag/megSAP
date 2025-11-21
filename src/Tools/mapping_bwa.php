<?php

/** 
	@page mapping_bwa
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_bwa", "Maps paired-end reads to a reference genome using bwa.");
$parser->addInfile("in1",  "Input file in FASTQ format. Forward read.", false);
$parser->addInfile("in2",  "Input file in FASTQ format. Reverse read.", false);
$parser->addOutfile("out",  "Output file in BAM format (sorted).", false);
//optional
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "GRCh38");
$parser->addString("sample", "Sample name to use in BAM header. If unset the basename of the 'out' file is used.", true, "");
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
$parser->addFlag("dedup", "Mark duplicates after alignment.");
$parser->addFlag("use_bwa1", "Use bwa instead of BWA-mem2.");
extract($parser->parse($argv));

//remove out file if existing
if(file_exists($out)) unlink($out);

//set read group information
if ($sample=="") $sample = basename2($out);
$group_props = array();
$group_props[] = "ID:{$sample}";
$group_props[] = "SM:{$sample}";
$group_props[] = "LB:{$sample}";
$group_props[] = "DT:".date("c");
$group_props[] = "PL:ILLUMINA";
if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $sample, false);
	if ($psample_info!=NULL)
	{
		$group_props[] = "PM:".$psample_info['device_type'];
		$group_props[] = "en:".$psample_info['sys_name'];
	}
}

//set input files for bwa-mem2 container
$genome = genome_fasta($build);

$in_files = [$in1, $in2, $genome];

//mapping with bwa
$pipeline = array();
$bwa_params = "mem $genome -K 100000000 -Y -R '@RG\\t".implode("\\t", $group_props)."' -t $threads -v 2";

//select the correct binary
if ($use_bwa1 || get_path("use_bwa1")) 
{
	$pipeline[] = array("", $parser->execApptainer("bwa", "bwa", "$bwa_params $in1 $in2", $in_files, [], true));
}
else
{
	$suffix = trim(get_path("bwa_mem2_suffix", false));
	$pipeline[] = ["", $parser->execApptainer("bwa-mem2", "bwa-mem2{$suffix}", "$bwa_params $in1 $in2", $in_files, [], true)];
}

//duplicate removal with samblaster
if ($dedup)
{
	$pipeline[] = ["", $parser->execApptainer("samblaster", "samblaster", "", [], [], true)];
}

//convert sam to bam with samtools
$tmp_unsorted = $parser->tempFile("_unsorted.bam");
$pipeline[] = ["", $parser->execApptainer("samtools", "samtools", "view -Sb - > $tmp_unsorted", [], [], true)];

//execute (BWA -> samblaster -> BAM conversion)
$parser->execPipeline($pipeline, "mapping");
$parser->sortBam($tmp_unsorted, $out, $threads, $build);
$parser->indexBam($out, $threads);

?>
