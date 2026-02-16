<?php

/** 
	@page analyze_dragen
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("analyze_dragen", "Performs mapping and variant calling using the Illumina DRAGEN server.");
$parser->addInfile("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);

//optional
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'name' is a valid processed sample name).", true);
$steps_all = array("ma", "vc", "cn", "sv", "re", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, cn=copy-number analysis, sv=structural-variant analysis, db=import into NGSD.", true, "ma,vc,cn,sv,re,db");
$parser->addInt("threads", "The maximum number of threads used for later megSAP analysis.", true, 2);
$parser->addString("rna_sample", "Processed sample name of the RNA sample which should be used for annotation.", true, "");
$parser->addFlag("no_queuing", "Do not queue megSAP analysis afterwards.");
$parser->addFlag("mapping_only", "Only map the data and remove variant calling.");
$parser->addFlag("debug", "Add debug output to the log file.");
$parser->addFlag("high_priority", "Queue megSAP analysis with high priority.");
$parser->addFlag("somatic", "Queue megSAP analysis with somatic flag.");
$parser->addString("user", "User used to queue megSAP analysis (has to be in the NGSD).", true, "");
$parser->addFlag("high_mem", "Run DRAGEN analysis in high memory mode for deep samples. (Also increases timeout to prevent job from being terminated.)");
$parser->addFlag("trim", "Perform adapter and quality trimming. Optional since trimming is usually done during the demultiplexing.");

extract($parser->parse($argv));

// ********************************* init *********************************//

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($folder."/analyze_dragen_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//determine processing system
$system_created_from_ngsd = (is_null($system) || $system=="");
$sys = load_system($system, $name);
$is_wes_or_panel = $sys['type']=="WES" || $sys['type']=="Panel" || $sys['type']=="Panel Haloplex";
$is_wgs = $sys['type']=="WGS";
$build = $sys['build'];
$genome = genome_fasta($build, false);
$dragen_path = "/opt/dragen/".get_path("dragen_version")."/bin/";

//stop on unsupported processing systems:
if (!$is_wes_or_panel && !$is_wgs) trigger_error("Can only analyze WGS/WES/panel samples not ".$sys['type']." samples!", E_USER_ERROR);

//sample checks:
if (in_array($sys['umi_type'], [ "MIPs", "ThruPLEX", "Safe-SeqS", "QIAseq", "IDT-xGen-Prism", "Twist"])) trigger_error("UMI handling is not supported by the DRAGEN pipeline! Please use local mapping instead.", E_USER_ERROR);
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
// explicitly set file privileges, since mkdir doesn't seem to work 
$parser->exec("chmod", "777 {$working_dir}");

//get ORA files
$files_forward = glob("{$folder}/*_L00[0-9]_R1_00[0-9].fastq.ora");
$files_reverse = glob("{$folder}/*_L00[0-9]_R2_00[0-9].fastq.ora");

//check number of ORAs for forward/reverse is equal
if (count($files_forward)!=count($files_reverse)) 
{
	trigger_error("Found mismatching forward and reverse read file count!\n Forward: {$folder}/*_L00[0-9]_R1_00[0-9].fastq.ora\n Reverse: {$folder}/*_L00[0-9]_R2_00[0-9].fastq.ora.", E_USER_ERROR);
}

//fallback to FASTQ
if (count($files_forward) == 0)
{
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

//fallback to BAM/CRAM
if (count($files_forward) == 0)
{
	$input_bam = $folder."/".$name.".cram";
	$input_index = $folder."/".$name.".cram.crai";
	if (!file_exists($input_bam))
	{
		$input_bam = $folder."/".$name.".bam";
		$input_index = $folder."/".$name.".bam.bai";
		if (!file_exists($input_bam))
		{
			trigger_error("Neither ORA, FASTQ, CRAM nor BAM file found in sample folder! Cannot perform analysis!", E_USER_ERROR);
		}
	}
	
	//convert BAM/CRAM to FASTQ
	//NOTE: it is not possible to override the sample ID when using BAM/CRAM as input, so external samples get the wrong sample IDs, which end up in the VCF files and cause errors in subsequent secondary analyses, e.g. trio analysis. The workround for this problem is to convert the data to FASTQ.
	trigger_error("Converting BAM/CRAM to FASTQ. This is slow and should not be done on the DRAGEN server. If possible, do the conversion before calling this script!", E_USER_NOTICE);
	$in_fq_for = $folder."/{$name}_BamToFastq_R1_001.fastq.gz";
	$in_fq_rev = $folder."/{$name}_BamToFastq_R2_001.fastq.gz";
	$tmp1 = $working_dir."/{$name}_BamToFastq_R1_001.fastq.gz";
	$tmp2 = $working_dir."/{$name}_BamToFastq_R2_001.fastq.gz";
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$input_bam} -ref {$genome} -out1 {$tmp1} -out2 {$tmp2}", [$genome, $folder, $working_dir]);
	$parser->moveFile($tmp1, $in_fq_for);
	$parser->moveFile($tmp2, $in_fq_rev);
	$files_forward = [$in_fq_for];
	$files_reverse = [$in_fq_rev];
}

//define output folder
$dragen_out_folder = $folder."/dragen/";

//check if input files are readable and output file is writeable
foreach ($files_forward as $in_file) 
{
	if (!is_readable($in_file)) trigger_error("Input file '{$in_file}' is not readable!", E_USER_ERROR);
}
foreach ($files_reverse as $in_file) 
{
	if (!is_readable($in_file)) trigger_error("Input file '{$in_file}' is not readable!", E_USER_ERROR);
}
$parser->exec("mkdir", "-p {$dragen_out_folder}");

// ********************************* call dragen *********************************//

$dragen_parameter = [];

//get library
$rglb = "UnknownLibrary";
$device_type = "UnknownDevice";
if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $name, false, true);
	$rglb = strtr($psample_info['sys_name'], ',', ' ');
	$device_type = $psample_info['device_type'];
}

//create FASTQ file list
$fastq_file_list_data = array("RGID,RGSM,RGLB,Lane,Read1File,Read2File,RGDT,RGPL");
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
	
	//extract lane
	for ($j=0; $j < count($filename_r1); $j++) 
	{ 
		if (starts_with("L0", $filename_r1[$j]))
		{
			if (starts_with("L0", $filename_r2[$j]))
			{
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

	$fastq_file_list_data[] = "{$name},{$name},{$rglb},{$lane},{$fastq_for},{$fastq_rev},".date("c").",{$device_type}";
}
$fastq_file_list = $parser->tempFile(".csv");
if ($debug) $parser->log("FastQ file list:", $fastq_file_list_data); //debug
file_put_contents($fastq_file_list, implode("\n", $fastq_file_list_data));
$dragen_parameter[] = "--fastq-list ".$fastq_file_list;

//create target region for exomes/panels
if ($is_wes_or_panel)
{
	$target_extended = $parser->tempFile("_roi_extended.bed");
	$target_region = file($sys['target_file'], FILE_IGNORE_NEW_LINES);
	$target_region_extended = array();
	foreach ($target_region as $line) 
	{
		if (starts_with($line, "#")) 
		{
			$target_region_extended[] = $line;
		}
		else
		{
			$parts = explode("\t", $line);
			if (count($parts) < 3) trigger_error("Invalid BED line: '{$line}'!", E_USER_ERROR);
			$parts[1] = intval($parts[1]) - 200;
			$parts[2] = intval($parts[2]) + 200;
			$target_region_extended[] = implode("\t", $parts);
		}	
	}
	file_put_contents($target_extended, implode("\n", $target_region_extended));
}

//trimming
if ($trim)
{
	//write adapters to fasta files
	$tmp_adapter_p5 = $parser->tempFile("_adapter_p5.fasta");
	$tmp_adapter_p7 = $parser->tempFile("_adapter_p7.fasta");
	file_put_contents($tmp_adapter_p5, implode("\n", [">header adapter1 p5", $sys["adapter1_p5"]]));
	file_put_contents($tmp_adapter_p7, implode("\n", [">header adapter2 p7", $sys["adapter2_p7"]]));

	$dragen_parameter[] = "--read-trimmers adapter,quality";
	$dragen_parameter[] = "--trim-adapter-read1 $tmp_adapter_p7";
	$dragen_parameter[] = "--trim-adapter-read2 $tmp_adapter_p5";
	$dragen_parameter[] = "--trim-min-quality 15";
}

//parameters
$dragen_parameter[] = "-r ".get_path("dragen_genome");
$dragen_parameter[] = "--ora-reference ".get_path("data_folder")."/dbs/oradata/";
$dragen_parameter[] = "--output-directory $working_dir";
$dragen_parameter[] = "--output-file-prefix {$name}";
$dragen_parameter[] = "--output-format CRAM"; //always use CRAM

//mapping:
$dragen_parameter[] = "--enable-map-align-output true";
$dragen_parameter[] = "--enable-cram-indexing true";
$dragen_parameter[] = "--enable-rh false"; //disabled RH special caller because this leads to variants with EVENTTYPE=GENE_CONVERSION that have no DP and AF entry and sometimes are duplicated (same variant twice in the VCF).
//set target region
if ($is_wes_or_panel) $dragen_parameter[] = "--vc-target-bed {$target_extended}";

if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $name, false, true);
}
//always mark duplicates
$dragen_parameter[] = "--enable-duplicate-marking true";

if (!$mapping_only)
{
	//small variant calling
	$dragen_parameter[] = "--enable-variant-caller true";
	$dragen_parameter[] = "--vc-min-base-qual 15"; //TODO remove when switching to DRAGEN 4.4
	//add gVCF output
	$dragen_parameter[] = "--vc-emit-ref-confidence GVCF";
	$dragen_parameter[] = "--vc-enable-vcf-output true";
	//disabled ML model because it leads to a sensitivity drop for Twist Exome V2 (see /mnt/storage2/users/ahsturm1/scripts/2025_03_21_megSAP_release_performance)
	if ($is_wes_or_panel)
	{
		$dragen_parameter[] = "--vc-ml-enable-recalibration false";
	}
	
	//CNVs
	if ($is_wgs)
	{
		$dragen_parameter[] = "--enable-cnv true";
		$dragen_parameter[] = "--cnv-enable-self-normalization true";
	}
	
	//SVs
	$dragen_parameter[] = "--enable-sv true";
	$dragen_parameter[] = "--sv-use-overlap-pair-evidence true"; //TODO remove when switching to DRAGEN 4.4
}

//high memory
if ($high_mem)
{
	$dragen_parameter[] = "--bin_memory 96636764160"; //90GB (default: 20 (GB), max on DRAGEN v4: 90)
	$dragen_parameter[] = "--vc-max-callable-region-memory-usage 41875931136"; //39GB (default: 13 (GB))
	$dragen_parameter[] = "--watchdog-active-timeout 3600"; //increase timeout to 1h (default: 10min)
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

//remove temporary created fastq files
foreach (glob("{$working_dir}/{$name}_BamToFastq_R?_00?.fastq.gz") as $tmp_fastq_file) 
{
	if ($debug) trigger_error("Removing tmp FASTQ file '{$tmp_fastq_file}'...", E_USER_NOTICE);
	unlink($tmp_fastq_file);
}

// ********************************* copy data to Sample folder *********************************//

$parser->log("Copying complete DRAGEN folder to output...");
//remove already existing output folder
if (file_exists($dragen_out_folder)) 
{
	$parser->exec("rm", "-rf $dragen_out_folder");
	mkdir($dragen_out_folder);
}
$parser->exec("cp", "-rf {$working_dir}/* {$dragen_out_folder}");

// ********************************* cleanup *********************************//

//delete working directory
$parser->exec("rm", "-rf $working_dir");

// ********************************* queue megSAP ****************************//


//parse megSAP parameter
$megSAP_args = array();
if (!$system_created_from_ngsd) $megSAP_args[] = "-system {$system}";
$megSAP_args[] = "-steps ".implode(",", $steps);
$megSAP_args[] = "-threads {$threads}";
if ($somatic) $megSAP_args[] = "-somatic";
if ($rna_sample != "") $megSAP_args[] = "-rna_sample {$rna_sample}";
$high_priority_str = ($high_priority)? "-high_priority " : ""; 

//queue analysis

if ($user == "") $user = "unknown";
$queuing_params = "-user {$user} -type 'single sample' -samples {$name} -ignore_running_jobs {$high_priority_str} -args '".implode(" ", $megSAP_args)."'";

if ($no_queuing) 
{
	trigger_error("megSAP queuing skipped! \nCommand to queue analysis: \n\tphp ".repository_basedir()."/src/Tools/db_queue_analysis.php {$queuing_params}", E_USER_NOTICE);
}
else 
{
	$parser->execTool("Tools/db_queue_analysis.php", $queuing_params);
}

//print to STDOUT executed successfully (because there is no exit code from SGE after a job has finished)
print "DRAGEN successfully finished!";

?>
