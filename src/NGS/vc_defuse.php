<?php
/**
 * @page fusion_detection_defuse
 * 
 * @todo remove SeqPurge
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// remove leading perl from defuse path
$tool_defuse = dirname(explode(" ",get_path("defuse"))[1]);

$parser = new ToolBase("vc_defuse", "Gene Fusion Detection using DeFuse");
$parser->addInfile("in_for", "FastQ file containing the forward reads", false, true);
$parser->addInfile("in_rev", "FastQ file containing the reverse reads", false, true);
$parser->addString("out", "Destination folder", false);
$parser->addString("sampleName", "Specify the sample name", false);
//optional parameters
$parser->addString("defuse_config", "Path to the defuse config file that should be used for the ananlysis", true, $tool_defuse."/config.txt");
$parser->addInt("numThreads", "Specifiy the number of threads used for the computations", true, 4);
$parser->addFlag("seqPurge", "Perform adapter and quality trimming using SeqPurge");
$parser->addInt("mq", "Minimal quality for quality trimming", true, 20);
$parser->addString("forAdapt", "Forward adapter for adapter clipping", true, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA");
$parser->addString("revAdapt", "Reverse adapter for adapter clipping", true, "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC");
$parser->addFlag("U", "Data is unzipped");

extract($parser->parse($argv));

$in1 = $in_for;
$in2 = $in_rev;

//run seqpurge first if specified
if($seqPurge) {
	
	if($U) { //if data is unzipped, we have to gzip it first since seqpurge only reads gzipped files
		$parser->log("Input data is unzipped, has to be gzipped first");
		$in1 = $parser->tempFile("_R1.fastq.gz");
		$in2 = $parser->tempFile("_R2.fastq.gz");
		$parser->log("Zipping input files for SeqPurge");
		$parser->exec("gzip", "-c ".$in_for." > ", true);
		$parser->exec("gzip", $in_rev, true);
	}
	
	//check quality encoding
	$parser->log("Checking quality encoding format of the input files");
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."FastqFormat", "-in $in_for", true);
	if (!contains($stdout[2], "Sanger"))
	{
		trigger_error("Input file '$in_for' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
	}
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."FastqFormat", "-in $in_rev", true);
	if (!contains($stdout[2], "Sanger"))
	{
		trigger_error("Input file '$in_rev' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
	}
	
	//run seq purge
	$parser->log("Performing adapter and quality trimming");
	$trim1 = $parser->tempFile("_R1_qt.bam");
	$trim2 = $parser->tempFile("_R2_qt.bam");
	$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 $in_for -in2 $in_rev -out1 $trim1 -out2 $trim2 -qcut $mq -a1 $forAdapt -a2 $revAdapt", true);
	
	//unzip trimmed fastq files for defuse <- only reads unzipped fastq files
	$parser->log("Data has to be unzipped after clipping and trimming to be applicable to defuse");
	$in1 = $parser->tempFile("_R1_qt_uz.fastq");
	$in2 = $parser->tempFile("_R2_qt_uz.fastq");
	$parser->exec("gunzip", "-c ".$trim1." > ".$in1, true);
	$parser->exec("gunzip", "-c ".$trim2." > ".$in2, true);
}

if(!$seqPurge && !$U) {
	//unzip fastq files for defuse
	$in1 = $parser->tempFile("_R1_uz.fastq");
	$in2 = $parser->tempFile("_R2_uz.fastq");
	$parser->exec("gunzip", "-c ".$in_for." > ".$in1, true);
	$parser->exec("gunzip", "-c ".$in_rev." > ".$in2, true);
}

$tempOutput = $parser->tempFolder("defuse");

//build command
$arguments = array();
$arguments[] = "-c $defuse_config";
$arguments[] = "-o $tempOutput";
$arguments[] = "-1 $in1";
$arguments[] = "-2 $in2";
$arguemnts[] = "-n $sampleName";
$arguments[] = "-p $numThreads";

//run defuse on adapter clipped and quality trimmed reads
$parser->exec(get_path("defuse"), implode(" ", $arguments), true);

//copy relevant files and delete the rest
$parser->exec("cp", $tempOutput."/*.tsv ".$out, true);
?>