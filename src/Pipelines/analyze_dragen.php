<?php

/** 
	@page analyze_dragen
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("analyze_dragen", "Maps paired-end reads to a reference genome and performs small/structural variant calling using the Illumina DRAGEN server.");
$parser->addInfile("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);

//optional
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'name' is a valid processed sample name).", true);
$steps_all = array("ma", "vc", "cn", "sv", "re", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, cn=copy-number analysis, sv=structural-variant analysis, db=import into NGSD.", true, "ma,vc,cn,sv,re,db");
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addString("rna_sample", "Processed sample name of the RNA sample which should be used for annotation.", true, "");
$parser->addFlag("queue_analysis", "Queue megSAP analysis afterwards.");
$parser->addFlag("enable_cnv", "Enable CNV calling with self-normalization (WGS only).");
$parser->addFlag("debug", "Add debug output to the log file.");

extract($parser->parse($argv));

// ********************************* init *********************************//

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($folder."/analyze_dragen_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//determine processing system
$sys = load_system($system, $name);
$is_wes = $sys['type']=="WES";
$is_wgs = $sys['type']=="WGS";
$is_panel = $sys['type']=="Panel" || $sys['type']=="Panel Haloplex";
$has_roi = $sys['target_file']!="";
$build = $sys['build'];
$genome = genome_fasta($build);

//check if valid reference genome is provided
$dragen_genome_path = get_path("dragen_genomes")."/".$build."/dragen/";
if (!file_exists($dragen_genome_path)) 
{
	trigger_error("Invalid genome build '".$build."' given. Path '".$dragen_genome_path."' not found on Dragen!", E_USER_ERROR);
}

//sample checks:
if (if (in_array($sys['umi_type'], [ "MIPs", "ThruPLEX", "Safe-SeqS", "QIAseq", "IDT-xGen-Prism", "Twist"]))) trigger_error("UMI handling is not supported by the DRAGEN pipeline! Please use local mapping instead.", E_USER_ERROR);
if ($sys['type']=="WGS (shallow)") trigger_error("Shallow genomes are not supported by the DRAGEN pipeline! Please use local mapping instead.", E_USER_ERROR);


//check for report config
if (db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
	list($rc_id) = report_config($db, $name);
	if($rc_id != -1) trigger_error("Report configuration for {$name} exists in NGSD! Cannot perform DRAGEN analysis!", E_USER_ERROR);
}

//remove mapping step:
if (in_array("ma", $steps))
{
	if (($key = array_search("ma", $steps)) !== false) unset($steps[$key]);
}
else
{
	trigger_error("Cannot perform DRAGEN analysis without mapping!", E_USER_ERROR);
}


//get ORA files

$files_forward = glob("{$folder}/*_L00[0-9]_R1_00[0-9].fastq.ora");
$files_reverse = glob("{$folder}/*_L00[0-9]_R2_00[0-9].fastq.ora");

//check number of ORAs for forward/reverse is equal
if (count($files_forward)!=count($files_reverse)) 
{
	trigger_error("Found mismatching forward and reverse read file count!\n Forward: {$folder}/*_L00[0-9]_R1_00[0-9].fastq.ora\n Reverse: {$folder}/*_L00[0-9]_R2_00[0-9].fastq.ora.", E_USER_ERROR);
}

if (count($files_forward) == 0)
{
	//fallback to FastQs
	//get FastQ files
	$files_forward = glob("{$folder}/*_R1_00[0-9].fastq.gz");
	$files_reverse = glob("{$folder}/*_R2_00[0-9].fastq.gz");
	

	//check for UMI index files
	$files_index = glob("{$folder}/*_index_*.fastq.gz");
	if (count($files_index) > 0) trigger_error("UMI handling is not supported by the DRAGEN pipeline! Please use local mapping instead.", E_USER_ERROR);

	//check number of FastQs for forward/reverse is equal
	if (count($files_forward)!=count($files_reverse)) 
	{
		trigger_error("Found mismatching forward and reverse read file count!\n Forward: {$folder}/*_R1_00[0-9].fastq.gz\n Reverse: {$folder}/*_R2_00[0-9].fastq.gz.", E_USER_ERROR);
	}	
}

if (count($files_forward) == 0)
{
	//fallback to BAM/CRAM
	$input_bam = $folder."/".$name.".cram";
	if (!file_exists($input_bam))
	{
		$input_bam = $folder."/".$name.".bam";
		if (!file_exists($input_bam))
		{
			trigger_error("Neither ORA, FastQ, CRAM nor BAM file found in sample folder! Cannot perform analysis!", E_USER_ERROR);
		}
	}
	//convert BAM/CRAM to FastQ
	$fastq_r1 = $parser->tempFile("{$name}_BamToFastq_R1_001.fastq.gz");
	$fastq_r2 = $parser->tempFile("{$name}_BamToFastq_R2_001.fastq.gz");
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$input_bam} -ref {$genome} -out1 {$fastq_r1} -out2 {$fastq_r2}", [$genome, $folder]);
	$files_forward = array($fastq_r1);
	$files_reverse = array($fastq_r2);

}


//define output files
$out_cram = $folder."/".$name.".cram";
$dragen_out_folder = $folder."/dragen_variant_calls/";
$parser->exec("mkdir", "-p {$dragen_out_folder}");
$out_vcf = $dragen_out_folder."/".$name."_dragen.vcf.gz";
$out_gvcf = $dragen_out_folder."/".$name."_dragen.gvcf.gz";
$out_sv = $dragen_out_folder."/".$name."_dragen_svs.vcf.gz";
$out_cnv = $dragen_out_folder."/".$name."_dragen_cnvs.vcf.gz";
$out_cnv_raw = $dragen_out_folder."/".$name."_dragen_cnvs.bw";




//check if input files are readable and output file is writeable
foreach ($files_forward as $in_file) 
{
	if (!is_readable($in1)) trigger_error("Input file '{$in_file}' is not readable!", E_USER_ERROR);
}
foreach ($files_reverse as $in_file) 
{
	if (!is_readable($in1)) trigger_error("Input file '{$in_file}' is not readable!", E_USER_ERROR);
}
if (!is_writable2($out_cram)) trigger_error("Output file '{$out_cram}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_vcf)) trigger_error("Output file '{$out_vcf}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_gvcf)) trigger_error("Output file '{$out_gvcf}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_sv)) trigger_error("Output file '{$out_sv}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_cnv)) trigger_error("Output file '{$out_cnv}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_cnv_raw)) trigger_error("Output file '{$out_cnv_raw}' is not writable!", E_USER_ERROR);


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

//create FastQ file list
$fastq_file_list_data = array("RGID,RGSM,RGLB,Lane,Read1File,Read2File");

//get library
$rglb = "UnknownLibrary";
if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $name, false, true);
	$rglb = $psample_info['sys_name'];
}

for ($i=0; $i < count($files_forward); $i++) 
{ 
	$fastq_for = $files_forward[$i];
	$fastq_rev = $files_reverse[$i];
	
	//extract lane
	$lane = 1;
	$filename_r1 = array_reverse(explode("_", basename2($fastq_for)));
	$filename_r2 = array_reverse(explode("_", basename2($fastq_rev)));

	//sanity check
	if (count($filename_r1)!=count($filename_r2)) trigger_error("Found mismatching file name in forward and reverse read file!\n Forward: ".implode("_", $filename_r1)."\n Reverse: ".implode("_", $filename_r2), E_USER_ERROR);
	
	for ($j=0; $j < count($filename_r1); $j++) 
	{ 
		if (starts_with("L0", $filename_r1[$j]))
		{
			if (starts_with("L0", $filename_r2[$j]))
			{
				//extract lane
				$lane_r1 = intval(substr($filename_r1[$j], 1));
				$lane_r2 = intval(substr($filename_r2[$j], 1));

				if($lane_r1 == $lane_r2)
				{
					$lane = $lane_r1;
					break;
				}
				else
				{
					trigger_error("Lane of forward and reverse read file do not match!\n Forward: ".implode("_", $filename_r1)."\n Reverse: ".implode("_", $filename_r2), E_USER_ERROR);
				}
			}
			else 
			{
				trigger_error("Found mismatching file name in forward and reverse read file!\n Forward: ".implode("_", $filename_r1)."\n Reverse: ".implode("_", $filename_r2), E_USER_ERROR);
			}
		}
	}

	//add line to CSV
	$fastq_file_list_data[] = "{$name},{$name},{$rglb},{$lane},{$fastq_for},{$fastq_rev}";
}
$fastq_file_list = $parser->tempFile(".csv");
if ($debug) $parser->log("FastQ file list:", $fastq_file_list_data); //debug
file_put_contents($fastq_file_list, implode("\n", $fastq_file_list_data));

//create target region for exomes/panels
if ($is_wes || $is_panel)
{
	$target_extended = $parser->tempFile("_roi_extended.bed");
	$parser->execApptainer("ngs-bits", "BedExtend", "-in ".$sys['target_file']." -n 200 -out $target_extended -fai {$genome}.fai", [$sys['target_file'], $genome]);

}

//parameters
$dragen_parameter = [];
$dragen_parameter[] = "-r ".$dragen_genome_path;
$dragen_parameter[] = "--fastq-list ".$fastq_file_list;
$dragen_parameter[] = "--output-directory $working_dir";
$dragen_parameter[] = "--output-file-prefix output";
$dragen_parameter[] = "--output-format CRAM"; //always use CRAM
$dragen_parameter[] = "--enable-map-align-output=true";
$dragen_parameter[] = "--enable-bam-indexing true";
$dragen_parameter[] = "--RGID $name";
$dragen_parameter[] = "--RGSM $name";
$dragen_parameter[] = "--RGDT ".date("c");
$dragen_parameter[] = "--enable-rh=false"; //disabled RH special caller because this leads to variants with EVENTTYPE=GENE_CONVERSION that have no DP and AF entry and sometimes are duplicated (same variant twice in the VCF).
if ($is_wgs)
{
	$dragen_parameter[] = "--enable-cnv true";
	$dragen_parameter[] = "--cnv-enable-self-normalization true";
	$dragen_parameter[] = "--vc-ml-enable-recalibration=false"; //disabled because it leads to a sensitivity drop for Twist Exome V2 (see /mnt/storage2/users/ahsturm1/scripts/2025_03_21_megSAP_release_performance)
}
//set target region
if ($is_wes || $is_panel) $dragen_parameter[] = "--vc-target-bed {$target_extended}";

if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $name, false, true);
	$dragen_parameter[] = "--RGPL '".$psample_info['device_type']."'";
	$dragen_parameter[] = "--RGLB '".$psample_info['sys_name']."'";
}
//always mark duplicates
$dragen_parameter[] = "--enable-duplicate-marking true";

//small variant calling
$dragen_parameter[] = "--enable-variant-caller true";
$dragen_parameter[] = "--vc-min-base-qual 15"; //TODO Marc re-validate with DRAGEN 4.4 (also for somatic)

//add gVCFs
$dragen_parameter[] = "--vc-emit-ref-confidence GVCF";
$dragen_parameter[] = "--vc-enable-vcf-output true";
//structural variant calling
$dragen_parameter[] = "--enable-sv true";
$dragen_parameter[] = "--sv-use-overlap-pair-evidence true"; //TODO Marc re-validate with DRAGEN 4.4
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

//copy CRAM/CRAI to sample folder
$parser->log("Copying CRAM/CRAI to output folder");
$parser->copyFile($working_dir."output.cram", $out);
$parser->copyFile($working_dir."output.cram.crai", $out.".crai");

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
if ($is_wgs)
{
	$parser->log("Copying CNVs to output folder");
	$parser->copyFile($working_dir."output.cnv.vcf.gz", $out_cnv);
	$parser->copyFile($working_dir."output.cnv.vcf.gz.tbi", $out_cnv.".tbi");
	$parser->copyFile($working_dir."output.target.counts.diploid.bw", $out_cnv_raw);
}

// ********************************* cleanup *********************************//

//delete working directory
$parser->exec("rm", "-rf $working_dir");

// ********************************* queue megSAP ****************************//

if ($queue_analysis)
{
	//parse megSAP parameter
	$megSAP_args = array();
	if ($system != "") $megSAP_args[] = "-system {$system}";
	$megSAP_args[] = "-steps ".implode(",", $steps);
	$megSAP_args[] = "-threads {$threads}";
	if ($rna_sample != "") $megSAP_args[] = "-rna_sample {$rna_sample}";

	//queue analysis
	$parser->execTool("Tools/queue_analysis.php", "-type 'single sample' -samples {$name} -args '".implode(" ", $megSAP_args)."'");
}

//print to STDOUT executed successfully (because there is no exit code from SGE after a job has finished)
print "DRAGEN successfully finished!";

?>
