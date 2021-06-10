<?php

/** 
	@page vc_manta_dragen
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_manta_dragen", "SV calling with Manta using the Illumina DRAGEN server.");
$parser->addInfile("in", "Input BAM file (with index).", false);
$parser->addOutfile("out", "Output VCF file containing structural variants.", false);
//optional
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "GRCh38");
$parser->addFlag("exome", "If set, manta settings for exome/panel analysis are used (no depth filtering).");
extract($parser->parse($argv));

// path to data folder and reference genomes
$dragen_data_path = get_path("dragen_data");
$dragen_genome_path = get_path("dragen_genomes");

//check if input files are readable and output file is writeable:
if (!is_readable($in)) trigger_error("Input file '.$in.' is not readable!", E_USER_ERROR);
if (!is_readable($in.".bai")) trigger_error("Input file '.$in.' has no index!", E_USER_ERROR);
if (!is_writable2($out)) trigger_error("Output file '.$out.' is not writable!", E_USER_ERROR);

// create temporary folder for analysis:
$working_dir = $dragen_data_path."megSAP_working_dir_SV/";

// delete folder from previous mapping
if (file_exists($working_dir)) $parser->exec("rm", "-rf $working_dir");

// create empty folder for analysis
if (!mkdir($working_dir, 0700)) trigger_error("Could not create working directory '".$working_dir."'!", E_USER_ERROR);

// copy BAM file to local storage
$local_bam = $working_dir.basename($in);
$local_bam_index = $working_dir.basename($in).".bai";
$parser->copyFile($in, $local_bam);
$parser->copyFile($in.".bai", $local_bam_index);

// generate parameter list for dragen
$dragen_parameter = array();
$dragen_parameter[] = "-r ".$dragen_genome_path.$build."/dragen/";
$dragen_parameter[] = "--bam-input ".$local_bam;
$dragen_parameter[] = "--enable-map-align false";
$dragen_parameter[] = "--enable-sv true";
if ($exome) $dragen_parameter[] = "--sv-exome true";
$dragen_parameter[] = "--output-directory $working_dir";
$dragen_parameter[] = "--output-file-prefix ".basename($out, ".vcf.gz");


// run SV calling
print basename($in)."\n";
$result = $parser->exec("dragen", implode(" ", $dragen_parameter));

// save stdout and stderr in output folder:
if ($result[2] != 0)
{
	// error during SV calling
	trigger_error("SV calling failed! \n   Stderr:\n".implode("\n", $result[1]), E_USER_ERROR);
}

// log mapping metrics:
$sv_metrics = file($working_dir.basename($out, ".vcf.gz").".sv_metrics.csv");
$parser->log("DRAGEN sv metrics:", $sv_metrics);

// fix rights
$parser->exec("chmod", "-R 755 $working_dir");


// copy result to sample folder
$parser->log("Copy alignment back to /mnt/...");
$parser->copyFile($working_dir.basename($out, ".vcf.gz").".sv.vcf.gz", $out);
$parser->copyFile($working_dir.basename($out, ".vcf.gz").".sv.vcf.gz.tbi", $out.".tbi");


// delete working directory
$parser->exec("rm", "-rf $working_dir");

// check if deletion was successful
if (file_exists($working_dir))
{
	trigger_error("Folder '$working_dir' could not be removed.", E_USER_ERROR);
}


?>
