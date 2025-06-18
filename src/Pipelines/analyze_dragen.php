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
$parser->addFlag("no_queuing", "Do not queue megSAP analysis afterwards.");
$parser->addFlag("dragen_only", "Perform only DRAGEN analysis and copy all output files without renaming (not compatible with later megSAP analysis).");
$parser->addFlag("debug", "Add debug output to the log file.");
$parser->addFlag("high_priority", "Queue megSAP analysis with high priority.");
$parser->addString("user", "User used to queue megSAP analysis (has to be in the NGSD).", true, "");

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
$is_wes = $sys['type']=="WES";
$is_wgs = $sys['type']=="WGS";
$is_panel = $sys['type']=="Panel" || $sys['type']=="Panel Haloplex";
$has_roi = $sys['target_file']!="";
$build = $sys['build'];
$genome = genome_fasta($build, false);

//stop on unsupported processing systems:
if (!($is_wes || $is_wgs || $is_panel)) trigger_error("Can only analyze WGS/WES/panel samples not ".$sys['type']." samples!", E_USER_ERROR);

//check if valid reference genome is provided
$dragen_genome_path = get_path("dragen_genomes")."/".$build."/dragen/";
if (!file_exists($dragen_genome_path)) 
{
	trigger_error("Invalid genome build '".$build."' given. Path '".$dragen_genome_path."' not found on Dragen!", E_USER_ERROR);
}

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

//disable megSAP queuing for DRAGEN only analyses
if ($dragen_only) $no_queuing = true;


if (!$no_queuing)
{
	//add vc, sv steps for later megSAP analysis
	if (!in_array("vc", $steps))
	{
		trigger_error("Small variant calling step missing in later megSAP analysis, adding step 'vc'!", E_USER_WARNING);
		$steps[] = "vc";
	}
	if (!in_array("sv", $steps))
	{
		trigger_error("Structural variant calling step missing in later megSAP analysis, adding step 'sv'!", E_USER_WARNING);
		$steps[] = "sv";
	}
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

$input_bam = "";
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
	//move BAM/CRAM to temp folder 
	if (file_exists("{$folder}/bams_for_mapping")) $parser->exec("rm", "-rf {$folder}/bams_for_mapping");
	$parser->exec("mkdir", "{$folder}/bams_for_mapping");
	$parser->moveFile($input_bam, "{$folder}/bams_for_mapping/".basename($input_bam));
	$input_bam = "{$folder}/bams_for_mapping/".basename($input_bam);

	/*
	//convert BAM/CRAM to FastQ (don't use temp since it is tiny on DRAGEN)
	$fastq_r1 = "{$working_dir}/{$name}_BamToFastq_R1_001.fastq.gz";
	$fastq_r2 = "{$working_dir}/{$name}_BamToFastq_R2_001.fastq.gz";
	//use samtools to avoid requirement of Apptainer
	$parser->exec("samtools", "fastq -1 {$fastq_r1} -2 {$fastq_r2} -0 /dev/null -s /dev/null -n {$input_bam} --reference {$genome} -@ 15");
	$files_forward = array($fastq_r1);
	$files_reverse = array($fastq_r2);
	*/
}

//define output files
$out_cram = $folder."/".$name.".cram";
$dragen_out_folder = $folder."/dragen_variant_calls/";
$out_vcf = $dragen_out_folder."/".$name."_dragen.vcf.gz";
$out_gvcf = $dragen_out_folder."/".$name."_dragen.gvcf.gz";
$out_sv = $dragen_out_folder."/".$name."_dragen_svs.vcf.gz";
$out_cnv = $dragen_out_folder."/".$name."_dragen_cnvs.vcf.gz";
$out_cnv_raw = $dragen_out_folder."/".$name."_dragen_cnvs.bw";

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
if (!is_writable2($out_cram)) trigger_error("Output file '{$out_cram}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_vcf)) trigger_error("Output file '{$out_vcf}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_gvcf)) trigger_error("Output file '{$out_gvcf}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_sv)) trigger_error("Output file '{$out_sv}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_cnv)) trigger_error("Output file '{$out_cnv}' is not writable!", E_USER_ERROR);
if (!is_writable2($out_cnv_raw)) trigger_error("Output file '{$out_cnv_raw}' is not writable!", E_USER_ERROR);


// ********************************* call dragen *********************************//

$dragen_parameter = [];

//get library
$rglb = "UnknownLibrary";
$device_type = "UnknownDevice";
if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $name, false, true);
	$rglb = $psample_info['sys_name'];
	$device_type = $psample_info['device_type'];
}

if ($input_bam != "")
{
	//use BAM/CRAM as input
	if (ends_with($input_bam, ".cram")) $dragen_parameter[] = "--cram-input {$input_bam}";
	else $dragen_parameter[] = "--bam-input {$input_bam}";

	// $dragen_parameter[] = "--RGID {$name}";
	// $dragen_parameter[] = "--RGSM {$name}";
	// $dragen_parameter[] = "--RGDT ".date("c");
	// $dragen_parameter[] = "--RGPL '{$rglb}'";
	// $dragen_parameter[] = "--RGLB '{$device_type}'";
}
else
{
	/*
	//concat Fastq files
	$fastq_r1 = "{$working_dir}/{$name}_merged_R1_001.fastq.gz";
	$fastq_r2 = "{$working_dir}/{$name}_merged_R2_001.fastq.gz";

	$parser->exec("cat", implode(" ", $files_forward)." > {$fastq_r1}");
	$parser->exec("cat", implode(" ", $files_reverse)." > {$fastq_r2}");
	*/
	
	//create FastQ file list
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

		//add line to CSV
		$fastq_file_list_data[] = "{$name},{$name},{$rglb},{$lane},{$fastq_for},{$fastq_rev},".date("c").",{$device_type}";

	}
	$fastq_file_list = $parser->tempFile(".csv");
	if ($debug) $parser->log("FastQ file list:", $fastq_file_list_data); //debug
	file_put_contents($fastq_file_list, implode("\n", $fastq_file_list_data));

	//add dragen parameter:
	$dragen_parameter[] = "--fastq-list ".$fastq_file_list;
	
	/*
	$dragen_parameter[] = "-1 ".$fastq_r1;
	$dragen_parameter[] = "-2 ".$fastq_r2;

	$dragen_parameter[] = "--RGID {$name}";
	$dragen_parameter[] = "--RGSM {$name}";
	$dragen_parameter[] = "--RGDT ".date("c");
	$dragen_parameter[] = "--RGPL '{$rglb}'";
	$dragen_parameter[] = "--RGLB '{$device_type}'";
	*/
}

//create target region for exomes/panels
if ($is_wes || $is_panel)
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

//parameters
$dragen_parameter[] = "-r ".$dragen_genome_path;
$dragen_parameter[] = "--ora-reference ".get_path("data_folder")."/dbs/oradata/";
$dragen_parameter[] = "--output-directory $working_dir";
$dragen_parameter[] = "--output-file-prefix {$name}";
$dragen_parameter[] = "--output-format CRAM"; //always use CRAM
$dragen_parameter[] = "--enable-map-align-output=true";
$dragen_parameter[] = "--enable-bam-indexing true";
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

//remove temporary created fastq files
foreach (glob("{$working_dir}/{$name}_BamToFastq_R?_00?.fastq.gz") as $tmp_fastq_file) 
{
	if ($debug) trigger_error("Removing tmp FASTQ file '{$tmp_fastq_file}'...", E_USER_NOTICE);
	unlink($tmp_fastq_file);
}

// ********************************* copy data to Sample folder *********************************//

if ($dragen_only)
{
	$parser->log("Copying complete DRAGEN folder to output...");
	//remove already existing output folder
	if (file_exists($dragen_out_folder)) 
	{
		$parser->exec("rm", "-rf $dragen_out_folder");
		mkdir($dragen_out_folder);
	}
	$parser->exec("cp", "-rf {$working_dir}/* {$dragen_out_folder}");
}
else
{
	//copy CRAM/CRAI to sample folder
	$parser->log("Copying CRAM/CRAI to output folder");
	$parser->copyFile($working_dir.$name.".cram", $out_cram);
	$parser->copyFile($working_dir.$name.".cram.crai", $out_cram.".crai");

	// copy small variant calls to sample folder
	$parser->log("Copying small variants to output folder");
	$parser->exec("mkdir", "-p {$dragen_out_folder}");
	$parser->copyFile($working_dir.$name.".hard-filtered.vcf.gz", $out_vcf);
	$parser->copyFile($working_dir.$name.".hard-filtered.vcf.gz.tbi", $out_vcf.".tbi");
	$parser->copyFile($working_dir.$name.".hard-filtered.gvcf.gz", $out_gvcf);
	$parser->copyFile($working_dir.$name.".hard-filtered.gvcf.gz.tbi", $out_gvcf.".tbi");

	// copy SV calls to sample folder
	$parser->log("Copying SVs to output folder");
	$parser->copyFile($working_dir.$name.".sv.vcf.gz", $out_sv);
	$parser->copyFile($working_dir.$name.".sv.vcf.gz.tbi", $out_sv.".tbi");

	// copy CNV calls
	if ($is_wgs)
	{
		$parser->log("Copying CNVs to output folder");
		$parser->copyFile($working_dir.$name.".cnv.vcf.gz", $out_cnv);
		$parser->copyFile($working_dir.$name.".cnv.vcf.gz.tbi", $out_cnv.".tbi");
		if (file_exists($working_dir.$name.".target.counts.diploid.bw")) $parser->copyFile($working_dir.$name.".target.counts.diploid.bw", $out_cnv_raw);
		else if (file_exists($working_dir.$name.".target.counts.bw")) $parser->copyFile($working_dir.$name.".target.counts.bw", $out_cnv_raw);
		else trigger_error("BigWig CNV file '{$working_dir}{$name}.target.counts.diploid.bw' not found!", E_USER_WARNING);
	}

}

// ********************************* cleanup *********************************//

//delete working directory
$parser->exec("rm", "-rf $working_dir");

// ********************************* queue megSAP ****************************//


//parse megSAP parameter
$megSAP_args = array();
if (!$system_created_from_ngsd) $megSAP_args[] = "-system {$system}";
$megSAP_args[] = "-steps ".implode(",", $steps);
$megSAP_args[] = "-threads {$threads}";
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
