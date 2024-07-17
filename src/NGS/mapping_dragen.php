<?php

/** 
	@page mapping_dragen
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_dragen", "Maps paired-end reads to a reference genome and performs small/structural variant calling using the Illumina DRAGEN server.");
$parser->addInfile("in1", "Input file in FASTQ format. Forward read.", false);
$parser->addInfile("in2", "Input file in FASTQ format. Reverse read.", false);
$parser->addOutfile("out", "Output BAM file (indexed).", false);
$parser->addOutfile("out_vcf", "Output VCF.GZ file for small variants (indexed).", false);
$parser->addOutfile("out_gvcf", "Output gVCF.GZ file for small variants (indexed).", false);
$parser->addOutfile("out_sv", "Output VCF.GZ file for structural variants (indexed).", false);
$parser->addOutfile("out_cnv", "Output VCF.GZ file for copy-number variants (indexed).", false);
$parser->addOutfile("out_cnv_raw", "Output BW file for copy-nubmer raw data.", false);
//optional
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "GRCh38");
$parser->addString("sample", "Sample name to use in BAM header. If unset the basename of the 'out' file is used.", true, "");
$parser->addFlag("dedup", "Mark duplicates after alignment.");
$parser->addFlag("enable_cnv", "Enable CNV calling with self-normalization (WGS only).");
$parser->addFlag("debug", "Add debug output to the log file.");
extract($parser->parse($argv));

// ********************************* init *********************************//

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
}

//if no sample name is given use output name
if ($sample=="") $sample = basename2($out);

//check if valid reference genome is provided
$dragen_genome_path = get_path("dragen_genomes")."/".$build."/dragen/";
if (!file_exists($dragen_genome_path)) 
{
	trigger_error("Invalid genome build '".$build."' given. Path '".$dragen_genome_path."' not found on Dragen!", E_USER_ERROR);
}

//check if input files are readable and output file is writeable
if (!is_readable($in1)) trigger_error("Input file '$in1' is not readable!", E_USER_ERROR);
if (!is_readable($in2)) trigger_error("Input file '$in2' is not readable!", E_USER_ERROR);
if (!is_writable2($out)) trigger_error("Output file '$out' is not writable!", E_USER_ERROR);
if (!is_writable2($out_vcf)) trigger_error("Output file '$out_vcf' is not writable!", E_USER_ERROR);
if (!is_writable2($out_sv)) trigger_error("Output file '$out_sv' is not writable!", E_USER_ERROR);

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
$dragen_parameter[] = "-r ".$dragen_genome_path;
$dragen_parameter[] = "-1 ".$in1;
$dragen_parameter[] = "-2 ".$in2;
$dragen_parameter[] = "--output-directory $working_dir";
$dragen_parameter[] = "--output-file-prefix output";
$dragen_parameter[] = "--enable-map-align-output=true";
$dragen_parameter[] = "--enable-bam-indexing true";
$dragen_parameter[] = "--RGID $sample";
$dragen_parameter[] = "--RGSM $sample";
$dragen_parameter[] = "--RGCN medical_genetics_tuebingen";
$dragen_parameter[] = "--RGDT ".date("c");
$dragen_parameter[] = "--vc-ml-enable-recalibration=false"; //disabled because it leads to a sensitivity drop for Twist Exome V2 SNVs of 0.5% (see /mnt/storage2/users/ahsturm1/scripts/2023_08_01_megSAP_performance/)
$dragen_parameter[] = "--enable-rh=false"; //disabled RH special caller because this leads to variants with EVENTTYPE=GENE_CONVERSION that have no DP and AF entry and sometimes are duplicated (same variant twice in the VCF).
if ($enable_cnv)
{
	$dragen_parameter[] = "--enable-cnv true";
	$dragen_parameter[] = "--cnv-enable-self-normalization true";
	$dragen_parameter[] = "--cnv-interval-width 1000";
}
if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $sample, false, true);
	$dragen_parameter[] = "--RGPL '".$psample_info['device_type']."'";
	$dragen_parameter[] = "--RGLB '".$psample_info['sys_name']."'";
}
if ($dedup) $dragen_parameter[] = "--enable-duplicate-marking true";

//small variant calling
$dragen_parameter[] = "--enable-variant-caller true";
$dragen_parameter[] = "--vc-min-read-qual 1";
$dragen_parameter[] = "--vc-min-base-qual 15";
//add gVCFs
$dragen_parameter[] = "--vc-emit-ref-confidence GVCF";
$dragen_parameter[] = "--vc-enable-vcf-output true";
//structural variant calling
$dragen_parameter[] = "--enable-sv true";
$dragen_parameter[] = "--sv-use-overlap-pair-evidence true";
$parser->log("DRAGEN parameters:", $dragen_parameter);

//run
$parser->exec("dragen_reset", "");
$parser->exec("LANG=en_US.UTF-8 dragen", implode(" ", $dragen_parameter)); //LANG is necessary to avoid the error "locale::facet::_S_create_c_locale name not valid" if the locale from the ssh source shell is not available on the Dragen server 
if ($debug)
{
	list($stdout) = exec2("ls $working_dir");
	$parser->log("working_dir content:", $stdout);
}

// ********************************* copy data back *********************************//

//copy BAM/BAI to sample folder
$parser->log("Copying BAM/BAI to output folder");
$parser->copyFile($working_dir."output.bam", $out);
$parser->copyFile($working_dir."output.bam.bai", $out.".bai");

// copy small variant calls to sample folder
$parser->log("Copying small variants to output folder");
$parser->copyFile($working_dir."output.hard-filtered.vcf.gz", $out_vcf);
$parser->copyFile($working_dir."output.hard-filtered.vcf.gz.tbi", $out_vcf.".tbi");
$parser->copyFile($working_dir."output.hard-filtered.gvcf.gz", $out_gvcf);
$parser->copyFile($working_dir."output.hard-filtered.gvcf.gz.tbi", $out_gvcf.".tbi");

// copy SV calls to sample folder
$parser->log("Copying SVs to output folder");
$parser->copyFile($working_dir."output.sv.vcf.gz", $out_sv);
$parser->copyFile($working_dir."output.sv.vcf.gz.tbi", $out_sv.".tbi");

// copy CNV calls
if ($enable_cnv)
{
	$parser->log("Copying CNVs to output folder");
	$parser->copyFile($working_dir."output.cnv.vcf.gz", $out_cnv);
	$parser->copyFile($working_dir."output.cnv.vcf.gz.tbi", $out_cnv.".tbi");
	$parser->copyFile($working_dir."output.target.counts.diploid.bw", $out_cnv_raw);
}

// ********************************* cleanup *********************************//

//delete working directory
$parser->exec("rm", "-rf $working_dir");

//print to STDOUT executed successfully (because there is no exit code from SGE after a job has finished)
print "DRAGEN successfully finished!";

?>
