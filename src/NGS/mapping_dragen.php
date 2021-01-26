<?php

/** 
	@page mapping_dragen
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_dragen", "Maps paired-end reads to a reference genome using the Illumina DRAGEN server.");
$parser->addInfile("in1", "Input file in FASTQ format. Forward read.", false);
$parser->addInfile("in2", "Input file in FASTQ format. Reverse read.", false);
$parser->addOutfile("out", "Output file in BAM format (sorted).", false);
//optional
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "GRCh37");
$parser->addString("sample", "Sample name to use in BAM header. If unset the basename of the 'out' file is used.", true, "");
$parser->addFlag("dedup", "Mark duplicates after alignment.");
extract($parser->parse($argv));


// check if valid reference genome is provided
if (!in_array($build, array("GRCh37", "GRCh38", "GRCh38_alt", "hg19")))
{
	trigger_error("Invalid genome build '".$build."' given!", E_USER_ERROR);
}

// if no sample name is given use output name
if ($sample=="") $sample = basename($out, ".bam");

// path to data folder and reference genomes
$dragen_data_path = get_path("dragen_data");
$dragen_genome_path = get_path("dragen_genomes");

//check if input files are readable and output file is writeable:
if (!is_readable($in1)) trigger_error("Input file '.$in1.' is not readable!", E_USER_ERROR);
if (!is_readable($in2)) trigger_error("Input file '.$in2.' is not readable!", E_USER_ERROR);
if (!is_writable2($out)) trigger_error("Output file '.$out.' is not writable!", E_USER_ERROR);

// create temporary folder for analysis:
$working_dir = $dragen_data_path."megSAP_working_dir/";

// delete folder from previous mapping
if (file_exists($working_dir)) $parser->exec("rm", "-rf $working_dir");

// create empty folder for analysis
if (!mkdir($working_dir, 0700)) trigger_error("Could not create working directory '".$working_dir."'!", E_USER_ERROR);

// copy FastQ files to local storage
$local_fastq1 = $working_dir.basename($in1);
$local_fastq2 = $working_dir.basename($in2);
$parser->copyFile($in1, $local_fastq1);
$parser->copyFile($in2, $local_fastq2);

// generate parameter list for dragen
$dragen_parameter = array();
$dragen_parameter[] = "-r ".$dragen_genome_path.$build."/dragen/";
$dragen_parameter[] = "-1 ".$local_fastq1;
$dragen_parameter[] = "-2 ".$local_fastq2;
$dragen_parameter[] = "--output-directory $working_dir";
$dragen_parameter[] = "--output-file-prefix ".basename($out, ".bam");
$dragen_parameter[] = "--enable-bam-indexing true";

// set read group information
$dragen_parameter[] = "--RGID $sample";
$dragen_parameter[] = "--RGSM $sample";
$dragen_parameter[] = "--RGCN medical_genetics_tuebingen";
$dragen_parameter[] = "--RGDT ".date("c");

if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $sample, false);
	$dragen_parameter[] = "--RGPL '".$psample_info['device_type']."'";
	$dragen_parameter[] = "--RGLB '".$psample_info['sys_name']."'";
}

// remove duplicates
if ($dedup) $dragen_parameter[] = "--enable-duplicate-marking true";


// run mapping
print "dragen ".implode(" ", $dragen_parameter)."\n";

$result = $parser->exec("dragen", implode(" ", $dragen_parameter));

// save stdout and stderr in output folder:
if ($result[2] != 0)
{
	// error during mapping
	trigger_error("Mapping failed! \n   Stderr:\n".implode("\n", $result[1]), E_USER_ERROR);
}

// log mapping metrics:
$mapping_metrics = file($working_dir.basename($out, ".bam").".mapping_metrics.csv");
$parser->log("DRAGEN mapping metrics:", $mapping_metrics);


// copy result to sample folder
$parser->log("Copy alignment back to /mnt/...");
$parser->copyFile($working_dir.basename($out), $out);
$parser->copyFile($working_dir.basename($out).".bai", $out.".bai");

// delete working directory
$parser->exec("rm", "-rf $working_dir");

// check if deletion was successful
if (file_exists($working_dir))
{
	trigger_error("Folder '$working_dir' could not be removed.", E_USER_ERROR);
}

print "\nDRAGEN mapping successfully finished, exit code 0.";
?>
