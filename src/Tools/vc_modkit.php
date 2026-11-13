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
$parser->addString("region", "Limit processing to specific region (e.g. for testing). Format: [CHR]:[START]-[END] or [CHR]", true);
$parser->addInfile("include_bed", "Limit processing to specific BED file", true);
$parser->addFlag("debug", "Enable modkit debug log.");

extract($parser->parse($argv));

//init
$genome = genome_fasta($build);
$log_file = $parser->tempFile("_modkit_pileup.log");
$tmp_folder = $parser->tempFolder("modkit");
//convert to absolute path
$bed = realpath2($bed);

//run modkit
$in_files = [];
$in_files[] = $bam;
$in_files[] = $genome;
$args = [];
$args[] = "pileup";
$args[] = $bam;
$args[] = $tmp_folder;
$args[] = "--phased";
$args[] = "--cpg";
$args[] = "--modified-bases 5mC 5hmC";
$args[] = "--combine-strands";
$args[] = "--ref ".$genome;
$args[] = "--threads ".$threads;
if ($debug) $args[] = "--log-filepath ".$log_file;
if (isset($region)) $args[] = "--region ".$region;
if (isset($include_bed)) 
{
	$args[] = "--include-bed ".realpath($include_bed);
	$in_files[] = realpath($include_bed);
}
$parser->execApptainer("modkit", "modkit", implode(" ", $args), $in_files);

//filter & copy logfile
if ($debug)
{
	$log = array_filter(file($log_file), fn($line) => !strpos($line, "[DEBUG]"));
	$parser->log("modkit pileup log file", $log);
}

//compress BED files
$bed_hp1 = dirname($bed)."/".basename($bed, ".bed.gz")."_hp1.bed.gz";
$bed_hp2 = dirname($bed)."/".basename($bed, ".bed.gz")."_hp2.bed.gz";
$parser->execApptainer("htslib", "bgzip", "-@ {$threads} -c -l 9 {$tmp_folder}/combined.bedmethyl > {$bed}", [], [dirname($bed)]);
$parser->execApptainer("htslib", "bgzip", "-@ {$threads} -c -l 9 {$tmp_folder}/hp1.bedmethyl > {$bed_hp1}", [], [dirname($bed_hp1)]);
$parser->execApptainer("htslib", "bgzip", "-@ {$threads} -c -l 9 {$tmp_folder}/hp2.bedmethyl > {$bed_hp2}", [], [dirname($bed_hp2)]);

//create index
// $parser->execApptainer("htslib", "tabix", "--help");
$parser->execApptainer("htslib", "tabix", "--verbosity 4 -p bed {$bed}", [$bed], [dirname($bed)], false, true, false, true);
$parser->execApptainer("htslib", "tabix", "--verbosity 4 -p bed {$bed_hp1}", [$bed_hp1], [dirname($bed_hp1)], false, true, false, true);
$parser->execApptainer("htslib", "tabix", "--verbosity 4 -p bed {$bed_hp2}", [$bed_hp2], [dirname($bed_hp2)], false, true, false, true);

//run summary
if (isset($summary))
{
	$log_file2 = $parser->tempFile("_modkit_summary.log");
	$args = [];
	$args[] = "summary";
	$args[] = $bam;
	$args[] = "--reference ".$genome;
	$args[] = "--threads ".$threads;
	if ($debug)	$args[] = "--log-filepath ".$log_file2;
	$args[] = " > ".$summary;
	$parser->execApptainer("modkit", "modkit", implode(" ", $args), $in_files, [dirname($summary)]);
	
	//debug logging
	if ($debug) $parser->log("modkit summary log file", file($log_file2));
}


?>
