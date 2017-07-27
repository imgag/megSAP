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
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "GRCh37");
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
$parser->addFlag("dedup", "Mark duplicates after alignment.");
extract($parser->parse($argv));

//remove out file if existing
if(file_exists($out)) unlink($out);

//set read group information
$group_props = array();
$basename = basename($out, ".bam");
$group_props[] = "ID:{$basename}";
$group_props[] = "SM:{$basename}";
$group_props[] = "LB:{$basename}";
$group_props[] = "CN:medical_genetics_tuebingen";
$group_props[] = "DT:".date("c");
$group_props[] = "PL:ILLUMINA";
$psample_info = get_processed_sample_info($basename,false);
if(!is_null($psample_info))
{
	$group_props[] = "PM:".$psample_info['device_type'];
	$group_props[] = "en:".$psample_info['sys_name'];
}

//mapping with bwa
$pipeline = array();
$bwa_params = "mem ".get_path("local_data")."/{$build}.fa -M -R '@RG\\t".implode("\\t", $group_props)."' -t $threads";
$pipeline[] = array(get_path("bwa"), "$bwa_params $in1 $in2");

//duplicate removal with samblaster
if ($dedup)
{
	$pipeline[] = array(get_path("samblaster"), "");
}

//convert sam to bam with samtools
$tmp_unsorted = $parser->tempFile();
$pipeline[] = array(get_path("samtools"), "view -1 - > $tmp_unsorted");

//execute (BWA -> samblaster -> BAM conversion)
$parser->execPipeline($pipeline, "mapping");

//sort BAM according to coordinates
$tmp_for_sorting = $parser->tempFile();
$parser->exec(get_path("samtools"), "sort -T $tmp_for_sorting -m 1G -@ ".min($threads, 4)." -o $out $tmp_unsorted", true);

//create index file
$parser->exec(get_path("samtools")." index", " $out", true);

?>
