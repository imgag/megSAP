<?php

/**
	@page mapping_star_htseq
	@todo move star files to temporary folder
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_star_htseq", "Alignment of RNASeq fastq files to a reference genome.");
$parser->addInfile("in1", "Input forward reads in FASTQ.GZ format.", false);
$parser->addOutfile("out", "Output BAM file name", false);

//optional arguments
$parser->addInfile("in2",  "Input Reverse reads in FASTQ.GZ format for paired-end alignment.", true);
$parser->addString("genome", "Path to dir with STAR reference genome index.", true, get_path("data_folder")."genomes/STAR/hg19");

$parser->addFlag("stranded", "Add this flag if library preparation was strand conserving.");
$parser->addFlag("U", "Fastq input files are uncompressed.");
$parser->addInt("p", "Number of parallel threads.", true, 4);

$parser->addString("llink", "Adapter forward read", true, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA");
$parser->addString("rlink", "Adapter reverse read", true, "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC");

$parser->addInt("Clip5p", "Remove n bases from 5' end of read.", true, 0);
$parser->addInt("Clip3p", "Remove n bases from 3' end of read.", true, 0);

$parser->addFlag("rmdup", "Remove PCR duplicates.");

$parser->addInt("Lmax", "Max length of read fraction for seed search", true, 50);

$parser->addFlag("noSplicing", "Prevent reads from getting spliced");
$parser->addFlag("wiggleOutput", "Output alignment in wiggle format.");

$parser->addFlag("longReads", "Use STAR version suitable for very long reads > 500nt.");

$parser->addFlag("useSharedMemory", "Load STAR Index into shared memory in order to run multiple STAR instances on the same genome.");
$parser->addFlag("readCounting", "Perform htseq-count read counting after mapping.");
$parser->addFlag("clipping", "Perform adapter clipping during mapping.");
$parser->addFlag("keepUnmapped", "Save unmapped reads as fasta files.");

extract($parser->parse($argv));

$outdir = realpath(dirname($out))."/";
$prefix = $outdir.basename($out, ".bam");

// Temporary prefix where STAR stores his intermediate files. The directory has to be removed before running STAR, otherwise STAR complains.
$STAR_tmp_folder = $parser->tempFolder();
$parser->exec("rm -r ", $STAR_tmp_folder, true);

//build command
$arguments = array();
$arguments[] = "--genomeDir $genome";
$arguments[] = "--readFilesIn $in1";

//for paired-end mapping
if(isset($in2))
{
  $arguments[] = "$in2";
}

$arguments[] = "--outFileNamePrefix $prefix";
$arguments[] = "--outTmpDir $STAR_tmp_folder";

$arguments[] = "--runThreadN $p";
$arguments[] = "--outBAMsortingThreadN $p";

if($clipping) {
	$arguments[] = "--clip3pAdapterSeq $llink $rlink";
}

$arguments[] = "--clip3pNbases $Clip3p --clip5pNbases $Clip5p";
$arguments[] = "--outSAMattributes All";

if(!$stranded) {
  $arguments[] = " --outSAMstrandField intronMotif";
}

if(!$U) {
  $arguments[] = "--readFilesCommand zcat";
}

// Options for chimeric alignment detection.                                                                                                                                                         
$arguments[] = "--chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax parameter 3";
$arguments[] = "--seedSearchStartLmax $Lmax";

// 2-pass mode is only possible if shared memory is not used.
if($useSharedMemory) {
	$arguments[] = "--genomeLoad LoadAndRemove";
} else {
	$arguments[] = "--genomeLoad NoSharedMemory";
	$arguments[] = "--twopassMode Basic";
}

if($noSplicing) {
	//min has to be larger than max, to prevent spliced mapping
	$arguments[] = "--alignIntronMin 2";
	$arguments[] = "--alignIntronMax 1";
}

$arguments[] = "--outSAMtype BAM SortedByCoordinate";
//give at least 20GB of RAM for sorting
$arguments[] = "--limitBAMsortRAM 21474836480";

//add a read group line to the final bam file
$arguments[] = "--outSAMattrRGline ID:1 PL:illumina PU:RGPU LB:N SM:RGSM CN:MFT";

// produces {$prefix}.Unmapped.out.mate1 and possibly {$prefix}.Unmapped.out.mate2
if($keepUnmapped) {
	$arguments[] = "--outReadsUnmapped Fastx";
} else {
	$arguments[] = "--outReadsUnmapped None";
}

if($wiggleOutput) {
	$arguments[] = "--outWigType wiggle";
}

if($readCounting) {
	$arguments[] = "--quantMode GeneCounts";
}

//run STAR
if($longReads)
{
	$parser->exec(get_path("STAR")."long", implode(" ", $arguments), true);
}
else
{
	$parser->exec(get_path("STAR"), implode(" ", $arguments), true);
}

// Write the final log file into the tool log and delete the final log file afterwards
$final_log = "{$prefix}Log.final.out";
$final_log_lines = file($final_log);
$parser->log("Final STAR log", $final_log_lines);

// Compress fasta files containing unmapped reads.
if (file_exists("{$prefix}Unmapped.out.mate1"))
{
	$parser->exec("gzip", "--suffix .fasta.gz {$prefix}Unmapped.out.mate1", true);
}
if (file_exists("{$prefix}Unmapped.out.mate2"))
{
	$parser->exec("gzip", "--suffix .fasta.gz {$prefix}Unmapped.out.mate2", true);
}

//output BAM file of STAR
$bam_mapped = "{$prefix}Aligned.sortedByCoord.out.bam";

//remove duplicates
if($rmdup)
{
	$parser->execTool("NGS/remove_duplicates.php", "-in $bam_mapped -out $bam_rmdup");
}
else
{
	rename($bam_mapped, $out);
}
$parser->exec(get_path("samtools"), "index $out", true);

//mapping QC
$stafile2 = "{$prefix}_stats_map.qcML";
$parser->exec(get_path("ngs-bits")."MappingQC", "-in $out -out $stafile2 -rna", true);

// Rename counts file if the read counting was enabled
if($readCounting)
{
	$parser->exec("mv", "{$prefix}ReadsPerGene.out.tab {$prefix}.HTSeq.tsv", false);
}

// We have to recreate the STAR temporary directory which was deleted by STAR. Otherwise the tool tests will complain that the directory does not exist and the test fails
mkdir($STAR_tmp_folder);

//cleanup (keep Chimeric.out.junction for downstream analysis)
$parser->exec("rm ", "{$prefix}Chimeric.out.sam", false);
$parser->exec("rm", "{$prefix}*Log.out", false);
$parser->exec("rm", "{$prefix}*Log.progress.out", false);
$parser->exec("rm", "{$prefix}SJ.out.tab", false);
$parser->exec("rm -r", "{$prefix}_STARgenome", false);
$parser->exec("rm -r", "{$prefix}_STARpass1", false);
$parser->exec("rm", $final_log, false);

?>
