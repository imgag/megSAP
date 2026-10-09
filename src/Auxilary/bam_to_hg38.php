<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("bam_to_hg38", "Converts a arbitrary BAM file to a BAM file for the current HG38 genome of megSAP.");
$parser->addInfile("in",  "Input BAM/CRAM file.", false);
$parser->addInfile("in_ref",  "Reference genome the input BAM/CRAM file is based on.", false);
$parser->addOutfile("out",  "Output BAM file name.", false);
$parser->addString("sample", "Sample name to use in BAM header. If unset, the read goup ID of the input BAM is used.", true, "");
$parser->addInt("threads", "Maximum number of threads used.", true, 6);
extract($parser->parse($argv));

//determine sample ID
if ($sample=="")
{
	list($header) = $parser->execApptainer("samtools", "samtools view", "-H $in", [$in], []);
	foreach($header as $line)
	{
		$line = trim($line);
		if (!starts_with($line, "@RG\t")) continue;
		
		$parts = explode("\t", $line);
		foreach($parts as $part)
		{
			if (starts_with($part, "ID:")) $sample = substr($part, 3);
		}
	}
}
if ($sample=="") trigger_error("Sample ID could not be determined from @RG header. Please provide it!", E_USER_ERROR);

//convert input BAM to FASTQs
$genome = genome_fasta("GRCh38");
$tmp_fq1 = $parser->tempFile("_1.fastq.gz");
$tmp_fq2 = $parser->tempFile("_2.fastq.gz");
$parser->execApptainer("ngs-bits", "BamToFastq", "-in $in -ref $in_ref -out1 $tmp_fq1 -out2 $tmp_fq2", [$in, $in_ref], []);

//re-map FASTQs to HG38
$out = realpath2($out);
$parser->execTool("Tools/mapping_bwa.php", "-in1 $tmp_fq1 -in2 $tmp_fq2 -dedup -build GRCh38 -out $out -threads $threads", [], [dirname($out)]);

?>