<?php

/** 
	@page vc_dragen_somatic
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_dragen_somatic", "Performs small/structural variant calling for tumor-normal samples using the Illumina DRAGEN server.");
$parser->addInfile("t_bam", "Input file - tumor sample BAM or CRAM.", false);
$parser->addInfile("n_bam", "Input file - normal sample BAM or CRAM.", false);
$parser->addString("prefix", "Prefix for output.", false);
$parser->addString("dragen_out_folder", "Folder for output.", false);

//optional
$parser->addInfile("normal_snvs", "Dragen SNV calls of the normalsample. Necessary if -out_cnv is given.", true);
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "GRCh38");

$parser->addFlag("is_targeted", "If the sequencing is a targeted sequencing.");
$parser->addFlag("high_mem", "");
$parser->addFlag("debug", "Add debug output to the log file.");

extract($parser->parse($argv));

//log server again as there can be multiple dragen servers
$parser->logServerEnvronment();

// ********************************* init *********************************//
$dragen_path = "/opt/dragen/".get_path("dragen_version")."/bin/";

if ($debug)
{
	list($stdout) = exec2("hostname");
	$parser->log("server:", $stdout);
	
	list($stdout) = exec2("whoami");
	$parser->log("user:", $stdout);
	
	list($stdout) = exec2("ulimit -a");
	$parser->log("limits:", $stdout);
	
	list($stdout) = exec2("env");
	$parser->log("env:", $stdout);

	list($stdout) = exec2("{$dragen_path}dragen_lic");
	$parser->log("dragen licence status before:", $stdout);
}

//check if input files are readable and output file is writeable
if (!is_readable($t_bam)) trigger_error("Input file '$t_bam' is not readable!", E_USER_ERROR);
if ($n_bam != "" && !is_readable($n_bam)) trigger_error("Input file '$n_bam' is not readable!", E_USER_ERROR);


// create empty folder for analysis
$working_dir = get_path("dragen_data")."/megSAP_working_dir/";
if (file_exists($working_dir))	
{
	$parser->exec("rm", "-rf $working_dir");
}
if (!mkdir($working_dir, 0777))
{
	trigger_error("Could not create working directory '".$working_dir."'!", E_USER_ERROR);
}

// ********************************* call dragen *********************************//

//parameters
$dragen_parameter = [];
$dragen_parameter[] = "-r ".get_path("dragen_genome");
if (ends_with($t_bam, ".bam"))
{
	$dragen_parameter[] = "--tumor-bam-input ".$t_bam;
}
else if (ends_with($t_bam, ".cram"))
{
	$dragen_parameter[] = "--tumor-cram-input ".$t_bam;
}
if ($n_bam != "")
{
	if (ends_with($n_bam, ".bam"))
	{
		$dragen_parameter[] = "--bam-input ".$n_bam;
	}
	else if (ends_with($n_bam, ".cram"))
	{
		$dragen_parameter[] = "--cram-input ".$n_bam;
	}
}
$dragen_parameter[] = "--output-directory $working_dir";
$dragen_parameter[] = "--output-file-prefix {$prefix}";
$dragen_parameter[] = "--enable-map-align false"; # cannot map multiple (tumor, normal) inputs at once
$dragen_parameter[] = "--pair-by-name true";

//parameters for  high memory mode
if ($high_mem)
{
	$dragen_parameter[] = "--bin_memory 85899345920"; //80GB (default: 20 (GB), max on DRAGEN v4: 90)
	$dragen_parameter[] = "--vc-max-callable-region-memory-usage 27917287424"; //26GB (default: 13 (GB))
}


//small variant calling
$dragen_parameter[] = "--enable-variant-caller true";
$dragen_parameter[] = "--vc-min-tumor-read-qual 3"; #default 3 for t-n, 20 for t-only 
$dragen_parameter[] = "--vc-min-base-qual 15";
$dragen_parameter[] = "--vc-callability-tumor-thresh 15"; # default 15 - minimum coverage in bam to try variant calling at that position
$dragen_parameter[] = "--vc-callability-normal-thresh 5"; # default 5
$dragen_parameter[] = "--vc-enable-unequal-ntd-errors false"; # disables model to correct FFPE errors.. //TODO get to work with model

//CNV calling:
// if ($normal_snvs == "") trigger_error("If outfile 'out_cnv' is given 'normal_snvs' have to be given as well.", E_USER_ERROR);

// $dragen_parameter[] = "--enable-cnv true";
// $dragen_parameter[] = "--cnv-somatic-enable-het-calling true";
// $dragen_parameter[] = "--cnv-normal-b-allele-vcf $normal_snvs";


//SV calling:
$dragen_parameter[] = "--enable-sv true";
$dragen_parameter[] = "--sv-use-overlap-pair-evidence true";
if ($is_targeted)
{
	$dragen_parameter[] = "--sv-exome true";
}

//MSI: microsattelite instability:
$msi_ref = get_path("data_folder") . "/dbs/DRAGEN/". ($is_targeted ? "WES" : "WGS") ."_v1.1.0_".$build."_microsatellites.list";
if (! is_file($msi_ref))
{
	trigger_error("Skipping dragen MSI calling because MSI reference: '{$msi_ref}' doesn't exist.", E_USER_WARNING);
}

if (is_file($msi_ref))
{
	$dragen_parameter[] = "--msi-command tumor-normal";
	$dragen_parameter[] = "--msi-microsatellites-file $msi_ref";
	$dragen_parameter[] = "--msi-coverage-threshold 60"; //recommended value is 60 for solid and 500 for liquid tumor (Dragen V4.2)
}

$parser->log("DRAGEN parameters:", $dragen_parameter);

//run
$parser->exec("{$dragen_path}dragen_reset", "");
$parser->exec("LANG=en_US.UTF-8 {$dragen_path}dragen", implode(" ", $dragen_parameter)); //LANG is necessary to avoid the error "locale::facet::_S_create_c_locale name not valid" if the locale from the ssh source shell is not available on the Dragen server 
if ($debug)
{
	list($stdout) = exec2("ls $working_dir");
	$parser->log("working_dir content:", $stdout);
}


// ********************************* copy data to Sample folder *********************************//

$parser->log("Copying complete DRAGEN folder to output...");
$parser->exec("cp", "-rf {$working_dir}/* {$dragen_out_folder}");



// ********************************* cleanup *********************************//

//delete working directory
$parser->exec("rm", "-rf $working_dir");

if ($debug)
{
	list($stdout) = exec2("{$dragen_path}dragen_lic");
	$parser->log("dragen licence status after:", $stdout);
}

//print to STDOUT executed successfully (because there is no exit code from SGE after a job has finished)
print "DRAGEN successfully finished!";

?>
