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
$parser->addOutfile("out",  "Output file in bam format (sorted).", false);
//optional
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "hg19");
$parser->addInfile("sheet",  "Illumina GAIIx sample sheet for read group generation.", true);
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
$parser->addFlag("dedup", "Marks duplicates after alignment.");
extract($parser->parse($argv));

//remove out file if existing
if(is_file($out)) unlink($out);

//determine read group information
$group_props = array();
if (!isset($sheet) || !file_exists($sheet)) 
{
	//from out (bam) file name, e.g.: GS120159_01.bam; nb: input files are currently temporary files without sample identifier.
	$basename = strtr(basename($out), array(".bam"=>""));
	list($gs, $ps) = explode("_", $basename."_99");
	$group_props[] = "ID:{$basename}";
	$group_props[] = "SM:{$gs}";
	$group_props[] = "LB:{$gs}_{$ps}";
}
else
{
	// from SampleSheet.csv:
	// FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
	// 64RKCAAXX,2,GS120159_01,hg19,TGACCA,-,N,2x78PE_MP,SJ,X-Chr
	$sheet_data = file($sheet);
	list($fcid, $lane, $ps, , $mid) = explode(",", $sheet_data[1]);
	$gs = $ps;
	if (strpos($ps, "_")!==false)
	{
		$gs = substr($ps, 0, strpos($ps, "_"));
	}
	$group_props[] = "ID:$ps.$fcid.$lane";
	$group_props[] = "SM:$gs";
	$group_props[] = "LB:$ps";
	$group_props[] = "PU:$fcid-$mid.$lane";
}
$group_props[] = "CN:medical_genetics_tuebingen";
$group_props[] = "DT:".date("c");
$group_props[] = "PL:ILLUMINA";

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
$parser->exec(get_path("ngs-bits")."BamIndex", "-in $out", true);

?>
