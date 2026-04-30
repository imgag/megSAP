<?php 
/** 
	@page vc_modkit
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_modkit", "ONT methylation annotation using modkit.");
$parser->addInfile("bam",  "Input file in BAM format. Note: .bam.bai file is required!", false);
$parser->addOutfile("bed", "Output BED file containing the methylation info. Bgzipped and indexed.", false);

//optional
$parser->addOutfile("summary", "Optional summary file.", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");

extract($parser->parse($argv));

//init
$genome = genome_fasta($build);
$log_file = $parser->tempFile("_modkit_pileup.log");
$tmp_folder = $parser->tempFolder("modkit");

//run modkit
$args = [];
$args[] = "pileup";
$args[] = $bam;
$args[] = $tmp_folder;
$args[] = "--phased";
$args[] = "--cpg";
$args[] = "--modified-bases 5mC 5hmC";
// $args[] = "--bgzf";
$args[] = "--ref ".$genome;
$args[] = "--log-filepath ".$log_file;
$args[] = "--threads ".$threads;

//set bind paths for modkit
$in_files = array();
$in_files[] = $bam;
$in_files[] = $genome;

$parser->execApptainer("modkit", "modkit", implode(" ", $args), $in_files);

//filter & copy logfile
$log = array_filter(file($log_file), fn($line) => !strpos($line, "[DEBUG]"));
$parser->log("modkit pileup log file", $log);

/*
//copy bed.gz
$parser->copyFile($tmp_folder."/combined.bed.gz", $bed);
$bed_hp1 = dirname($bed)."/".basename($bed, ".bed.gz")."_hp1.bed.gz";
$bed_hp2 = dirname($bed)."/".basename($bed, ".bed.gz")."_hp2.bed.gz";
$parser->copyFile($tmp_folder."/hp1.bed.gz", $bed_hp1);
$parser->copyFile($tmp_folder."/hp2.bed.gz", $bed_hp2);
*/

//compress BED files
$bed_hp1 = dirname($bed)."/".basename($bed, ".bed.gz")."_hp1.bed.gz";
$bed_hp2 = dirname($bed)."/".basename($bed, ".bed.gz")."_hp2.bed.gz";
$parser->execApptainer("htslib", "bgzip", "-@ {$threads} -c -l 9 {$tmp_folder}/combined.bedmethyl > {$bed}", [], [dirname($bed)]);
$parser->execApptainer("htslib", "bgzip", "-@ {$threads} -c -l 9 {$tmp_folder}/hp1.bedmethyl > {$bed_hp1}", [], [dirname($bed_hp1)]);
$parser->execApptainer("htslib", "bgzip", "-@ {$threads} -c -l 9 {$tmp_folder}/hp2.bedmethyl > {$bed_hp2}", [], [dirname($bed_hp2)]);

//create index
$parser->execApptainer("htslib", "tabix", "-f -p bed {$bed}", [$bed], [dirname($bed)]);
$parser->execApptainer("htslib", "tabix", "-f -p bed {$bed_hp1}", [$bed_hp1], [dirname($bed_hp1)]);
$parser->execApptainer("htslib", "tabix", "-f -p bed {$bed_hp2}", [$bed_hp2], [dirname($bed_hp2)]);

//run summary
if (isset($summary))
{
	$log_file2 = $parser->tempFile("_modkit_summary.log");
	$args = [];
	$args[] = "summary";
	$args[] = $bam;
	$args[] = "--threads ".$threads;
	$args[] = "--log-filepath ".$log_file2;
	$args[] = " > ".$summary;

	$parser->execApptainer("modkit", "modkit", implode(" ", $args), $in_files, [dirname($summary)]);
	//copy logfile 
	$parser->log("modkit summary log file", file($log_file2));
}


?>
