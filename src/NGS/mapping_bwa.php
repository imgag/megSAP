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
extract($parser->parse($argv));

//remove out file if existing
if(file_exists($out)) unlink($out);

//set read group information
if ($sample=="") $sample = basename2($out);
$group_props = array();
$group_props[] = "ID:{$sample}";
$group_props[] = "SM:{$sample}";
$group_props[] = "LB:{$sample}";
$group_props[] = "CN:medical_genetics_tuebingen";
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

//mapping with bwa
$pipeline = array();
$bwa_params = "mem ".genome_fasta($build, true, false)." -K 100000000 -Y -R '@RG\\t".implode("\\t", $group_props)."' -t $threads -v 2";

//select the correct binary
if (get_path("use_bwa1")) 
{
	$bwa_binary = get_path("bwa");
}
else
{
	$bwa_binary = get_path("bwa-mem2"); 
}
$pipeline[] = array($bwa_binary, "$bwa_params $in1 $in2");

//duplicate removal with samblaster
if ($dedup)
{
	$pipeline[] = array(get_path("samblaster"), "");
}

//convert sam to bam with samtools
$tmp_unsorted = $parser->tempFile("_unsorted.bam");
$pipeline[] = array(get_path("samtools"), "view -Sb - > $tmp_unsorted");

//execute (BWA -> samblaster -> BAM conversion)
$parser->execPipeline($pipeline, "mapping");
$parser->sortBam($tmp_unsorted, $out, $threads);
$parser->indexBam($out, $threads);

?>
