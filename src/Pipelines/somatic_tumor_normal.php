<?php

/**
 * @page somatic_tumor_normal
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_tumor_normal", "Differential analysis of tumor and normal DNA samples.");

//mandatory
$parser->addInfile("t_bam", "Tumor sample BAM file.", false);
$parser->addString("out_folder", "Output folder.", false);
$parser->addInfile("n_bam", "Normal sample BAM file.", false);

//optional
$parser->addInfile("t_rna_bam", "Tumor RNA sample BAM file.", true);
$parser->addString("prefix", "Output file prefix.", true, "somatic");

$steps_all = array("vc", "vi", "cn", "an", "msi", "an_rna", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\n" .
	"vc=variant calling (small variants and SVs), an=annotation (small variants and SVs),\n" .
	"cn=copy-number analysis, msi=microsatellite analysis,\n".
	"an_rna=annotate data from somatic RNA files,\n".
	"vi=virus detection, db=database import",
	true, "vc,cn,an,msi,vi,an_rna,db");

$parser->addInfile("system", "Processing system file used for tumor DNA sample (resolved from NGSD via tumor BAM by default).", true);
$parser->addInfile("n_system", "Processing system file used for normal DNA sample (resolved from NGSD via normal BAM by default).", true);
$parser->addFlag("skip_contamination_check", "Skips check of female tumor sample for male SRY DNA.");
$parser->addFlag("skip_correlation", "Skip sample correlation check.");
$parser->addFlag("skip_low_cov", "Skip low coverage statistics.");
$parser->addFlag("skip_signatures", "Skip calculation of mutational signatures.");
$parser->addFlag("skip_HRD", "Skip calculation HRD.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
$parser->addFlag("use_dragen", "Use Illumina dragen for somatic variant calling.");
$parser->addFlag("use_deepsomatic_test", "Use DeepSomatic for somatic variant calling. (normally set in settings.ini)");
$parser->addFlag("validation", "Option used for analyzing validation samples. Ignores checks: flagging in NGSD, report config, correlation");
//default cut-offs
$parser->addFloat("min_correlation", "Minimum correlation for tumor/normal pair.", true, 0.8);
$parser->addFloat("min_depth_t", "Tumor sample coverage cut-off for low coverage statistics.", true, 60);
$parser->addFloat("min_depth_n", "Normal sample coverage cut-off for low coverage statistics.", true, 60);
$parser->addInt("min_cov_files", "Minimum number of required tumor-normal pairs for CNV calling.", true, 7);
$parser->addString("cnv_baseline_pos", "baseline region for ClinCNV, format e.g. chr1:12-12532",true);
$parser->addFloat("cnv_wgs_purity_step", "ClinCNV step size while checking purity of clones",true);
$parser->addString("rna_ref_tissue", "Reference data for RNA annotation", true);
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));

###################################### AUXILARY FUNCTIONS ######################################

//Checks $baf_folder for missing B-AF files and creates them if neccessary
function complement_baf_folder($t_id, $n_id, $t_n_id_file, $baf_folder, &$db_conn, $build)
{
	global $parser;
	
	$ids = file($t_n_id_file);
	
	$already_in_list = false;
	//Check whether tumor-normal pair is already included in csv list in tumor coverage folder or not
	foreach($ids as $line)
	{
		if(starts_with($line,'#')) continue;
		if(empty(trim($line))) continue;
		list($t,$n) = explode(",",trim($line));
		if($t == $t_id && $n == $n_id)
		{
			$already_in_list = true;
			break;
		}
	}
	
	if (!$already_in_list) $ids[] = "{$t_id},{$n_id}";
	
	foreach($ids as $line)
	{
		if(starts_with($line,"#")) continue;
		list($tid,$nid) = explode(",",trim($line));
		if(!file_exists("{$baf_folder}/{$nid}.tsv"))
		{
			$ninfo = get_processed_sample_info($db_conn,$nid);
			$n_gsvar = $ninfo["ps_folder"] ."/{$nid}.GSvar";
			$n_bam = $ninfo["ps_bam"];
			$parser->execTool("Auxilary/create_baf_file.php", "-gsvar $n_gsvar -bam $n_bam -build $build -out_file {$baf_folder}/{$nid}.tsv");
		}
		if(!file_exists("{$baf_folder}/{$tid}.tsv"))
		{
			$ninfo = get_processed_sample_info($db_conn,$nid);
			$n_gsvar = $ninfo["ps_folder"] ."/{$nid}.GSvar";
			$tinfo = get_processed_sample_info($db_conn,$tid);
			$t_bam = $tinfo["ps_bam"];
			$parser->execTool("Auxilary/create_baf_file.php", "-gsvar $n_gsvar -bam $t_bam -build $build -out_file {$baf_folder}/{$tid}.tsv");
		}
	}
	
	//add it to the file after bafs were succesfully created if it wasn't in the list yet
	if(!$already_in_list)
	{
		file_put_contents($t_n_id_file,"{$t_id},{$n_id}\n", FILE_APPEND | LOCK_EX);
	}
}

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all))
	{
		trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	}
}

//check dragen requirements
if (in_array("vc", $steps)  && $use_dragen)
{
	$dragen_input_folder = get_path("dragen_in");
	$dragen_output_folder = get_path("dragen_out");
	if (!file_exists($dragen_input_folder))	
	{
		trigger_error("DRAGEN input folder \"".$dragen_input_folder."\" does not exist!", E_USER_ERROR);
	}
	if (!file_exists($dragen_output_folder))
	{
		trigger_error("DRAGEN input folder \"".$dragen_output_folder."\" does not exist!", E_USER_ERROR);
	}
}


###################################### SCRIPT START ######################################
//check which caller to use
$use_deepsomatic = $use_deepsomatic_test ?: get_path("use_deepsomatic");

if($validation)
{
	$skip_correlation = true;
	$skip_contamination_check = true;
	$skip_signatures = true;
}

if (!file_exists($out_folder))
{
	exec2("mkdir -p $out_folder");
}

$out_folder = realpath($out_folder);

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($out_folder."/somatic_tumor_normal_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//output prefix
$full_prefix = "{$out_folder}/{$prefix}";

//IDs, system and target region
$t_id = basename2($t_bam);
$t_basename = dirname($t_bam)."/".$t_id;
$sys = load_system($system, $t_id);
$roi = trim($sys["target_file"]);
$ref_genome = genome_fasta($sys['build']);

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);
}

//check that ROI is sorted
if ($roi!="")
{
	$roi = realpath($roi);
	if (!bed_is_sorted($roi)) trigger_error("Target region file not sorted: ".$roi, E_USER_ERROR);
}

//normal sample data
$n_id = basename2($n_bam);
$n_folder = dirname($n_bam);
$n_basename = dirname($n_bam)."/".$n_id;
$n_sys = load_system($n_system, $n_id);
$ref_folder_n = get_path("data_folder")."/coverage/".$n_sys['name_short'];

check_genome_build($t_bam, $sys['build']);
check_genome_build($n_bam, $n_sys['build']);

if ($sys["name_short"] != $n_sys["name_short"] && in_array("cn", $steps))
{
	trigger_error("Tumor and normal sample were sequenced with different processing systems - CNVs cannot be calculated. Removing 'cn' step!",E_USER_WARNING);
	
	$key = array_search("cn", $steps);
	unset($steps[$key]);
}

//Check whether both samples have same processing system
if($roi != $n_sys["target_file"])
{
	#test that tumor target is a subset of normal target
	exec($parser->execApptainer("ngs-bits", "BedSubtract", "-in ".realpath($roi)." -in2 ".realpath($n_sys["target_file"]), [$roi, $n_sys["target_file"]], [], true), $output, $return_var);
	
	foreach ($output as $line)
	{
		if ($line == "" || starts_with($line, "#")) continue;
		
		trigger_error("Tumor sample $t_id  target region is different from, and not a subset of, the normal sample $n_id target region.",E_USER_ERROR);
	}
	
	trigger_error("Tumor sample $t_id and normal sample $n_id have different target regions.",E_USER_WARNING);
}

//Abort if calling is requested and somatic report config exists in NGSD
if (db_is_enabled("NGSD") && !$validation)
{
	$db = DB::getInstance("NGSD", false);
	list($config_id, $config_vars_exist, $config_cnvs_exist, $config_svs_exists) = somatic_report_config($db, $t_id, $n_id);
	if (in_array("vc", $steps) && $config_vars_exist)
	{
		trigger_error("Somatic report configuration with SNVs exists in NGSD! Delete somatic report configuration for reanalysis of step 'vc'.", E_USER_ERROR);
	}
	if (in_array("vc", $steps) && $config_svs_exists)
	{
		trigger_error("Somatic report configuration with SVs exists in NGSD! Delete somatic report configuration for reanalysis of step 'vc'.", E_USER_ERROR);
	}
	if (in_array("cn", $steps) && $config_cnvs_exist)
	{
		trigger_error("Somatic report configuration with CNVs exists in NGSD! Delete somatic report configuration for reanalysis of step 'cn'.", E_USER_ERROR);
	}
}

//sample similarity check
$bams = array_filter([$t_bam, $n_bam]);

if ($skip_correlation || $validation)
{
	trigger_error("Genotype correlation check has been disabled!", E_USER_NOTICE);
}
else
{
	$in_files = $bams;
	$in_files[] = $ref_genome;
	$args_similarity = [
		"-in ".implode(" ", $bams),
		"-mode bam",
		"-build ".ngsbits_build($sys['build']),
		"-ref {$ref_genome}"
	];
	if (!empty($roi))
	{
		$args_similarity[] = "-roi ".realpath($roi);
		$in_files[] = $roi;
	}
	$output = $parser->execApptainer("ngs-bits", "SampleSimilarity", implode(" ", $args_similarity), $in_files);

	//extract colum 3 from output
	$table = array_map(
		function($str) { return explode("\t", $str); },
		array_slice($output[0], 1)
	);
	$correlation = array_column($table, 3);
	if (min($correlation) < $min_correlation)
	{
		trigger_error("Genotype correlation lower than {$min_correlation}!\n" . implode("\n", $output[0]), E_USER_ERROR);
	}
}

//check SRY coverage of tumor for female samples, can be a hint for contamination with male DNA

if (!$skip_contamination_check)
{
	$parser->execTool("Auxilary/check_sry_coverage.php", "-t_bam $t_bam -build ".$sys['build']);
}
else trigger_error("Skipping check of female tumor sample $t_bam for contamination with male genomic DNA.", E_USER_WARNING);

// Check samples are flagged correctly in NGSD
if( db_is_enabled("NGSD") && !$validation)
{
	$db = DB::getInstance("NGSD");
	
	$tinfo = get_processed_sample_info($db, $t_id, false);

	if(!is_null($tinfo) && $tinfo["is_tumor"] != 1)
	{
		trigger_error("Please check tumor processed sample {$t_id} in NGSD. The sample is not flagged as tumor tissue.", E_USER_ERROR);
	}
	
	$ninfo = get_processed_sample_info($db, $n_id, false);
	if(!is_null($ninfo) && $ninfo["is_tumor"] != 0)
	{
		trigger_error("Please check normal processed sample {$n_id} in NGSD. The sample is flagged as tumor tissue.", E_USER_ERROR);
	}
}

//low coverage statistics
$low_cov = "{$full_prefix}_stat_lowcov.bed";					// low coverage BED file
if ($sys['type'] !== "WGS" && !empty($roi) && !$skip_low_cov)
{
	$parser->execApptainer("ngs-bits", "BedLowCoverage", "-in ".realpath($roi)." -bam $t_bam -out $low_cov -cutoff $min_depth_t -threads {$threads} -ref {$ref_genome}", [$roi, $t_bam, $ref_genome], [dirname($low_cov)]);
	//combined tumor and normal low coverage files
	//normal coverage is calculated only for tumor target region
	$low_cov_n = $parser->tempFile("_nlowcov.bed");
	$parser->execApptainer("ngs-bits", "BedLowCoverage", "-in ".realpath($roi)." -bam $n_bam -out $low_cov_n -cutoff $min_depth_n -threads {$threads} -ref {$ref_genome}", [$roi, $n_bam, $ref_genome]);
	$parser->execPipeline([
		["", $parser->execApptainer("ngs-bits", "BedAdd", "-in $low_cov $low_cov_n", [$low_cov], [], true)],
		["", $parser->execApptainer("ngs-bits", "BedMerge", "-out $low_cov", [], [dirname($low_cov)], true)],
	], "merge low coverage BED files");
	// annotate with gene names
	if (db_is_enabled("NGSD"))
	{
		$parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-in $low_cov -extend 25 -out $low_cov", [$low_cov]);
	}
}

//variant calling
$manta_indels  = $full_prefix . "_manta_var_smallIndels.vcf.gz";	// small indels from manta
$manta_sv      = $full_prefix . "_manta_var_structural.vcf.gz";		// structural variants (vcf)
$manta_sv_bedpe= $full_prefix . "_manta_var_structural.bedpe"; 		// structural variants (bedpe)
$variants      = $full_prefix . "_var.vcf.gz";					// variants
$ballele       = $full_prefix . "_bafs.igv";					// B-allele frequencies
$hla_file_tumor = "{$t_basename}_hla_genotyper.tsv";
$hla_file_normal = "{$n_basename}_hla_genotyper.tsv";

if (in_array("vc", $steps))
{
	//get HLA types
	$parser->execTool("Tools/hla_genotyper.php", "-bam $t_bam -name $t_id -out ".$hla_file_tumor);
	$parser->execTool("Tools/hla_genotyper.php", "-bam $n_bam -name $n_id -out " . $hla_file_normal);
	
	// structural variant calling
	if (!$sys['shotgun'])
	{
		trigger_error("Structural variant calling deactivated for amplicon samples.", E_USER_NOTICE);
	}
	else if ($sys['umi_type'] === "ThruPLEX")
	{
		trigger_error("Structural variant calling deactivated for ThruPLEX samples.", E_USER_NOTICE);
	}
	else
	{
		$args_manta = [
			"-t_bam {$t_bam}",
			"-out {$manta_sv}",
			"-build ".$sys['build'],
			"-smallIndels {$manta_indels}",
			"-threads {$threads}"
		];
		$args_manta[] = "-bam $n_bam";

		if ($sys['type'] !== "WGS") //use exome flag for non targeted / exome samples (i.e. non WGS samples)
		{
			$args_manta[] = "-exome";
		}
		if (!empty($roi))
		{
			$args_manta[] = "-target {$roi}";
		}
		$parser->execTool("Tools/vc_manta.php", implode(" ", $args_manta));
		
		$parser->execApptainer("ngs-bits", "VcfToBedpe", "-in $manta_sv -out $manta_sv_bedpe", [$manta_sv], [dirname($manta_sv_bedpe)]);

		$parser->execTool("Tools/bedpe2somatic.php", "-in $manta_sv_bedpe -out $manta_sv_bedpe -tid $t_id -nid $n_id");
		
		if( db_is_enabled("NGSD") )
		{
			$parser->execApptainer("ngs-bits", "BedpeGeneAnnotation", "-in $manta_sv_bedpe -out $manta_sv_bedpe -add_simple_gene_names", [$manta_sv_bedpe]);
		}
	}
	
	if ($use_dragen)
	{
		list($server) = exec2("hostname -f");
		//DRAGEN OUTFILES
		$dragen_output_vcf = "$dragen_output_folder/{$prefix}_dragen.vcf.gz";
		$dragen_output_msi = "$dragen_output_folder/{$prefix}_dragen_msi.json";
		$dragen_output_svs = "$dragen_output_folder/{$prefix}_dragen_svs.vcf.gz";
		$dragen_output_cnvs = "$dragen_output_folder/{$prefix}_dragen_cnvs.vcf.gz";
		$dragen_log_file = "$dragen_output_folder/{$prefix}_dragen.log";
		$sge_logfile = date("YmdHis")."_".implode("_", $server)."_".getmypid();
		$sge_update_interval = 300; //5min
		
		$t_bam_dragen = $t_bam;
		$n_bam_dragen = $n_bam;
		if (contains($t_bam, "/tmp/"))
		{
			#running based on a local bam file -> copy it to be reachable from dragen:
			$t_bam_dragen = $dragen_input_folder.basename($t_bam);
			$parser->copyFile($t_bam, $t_bam_dragen);
		}
		
		if (contains($n_bam, "/tmp/"))
		{
			#running based on a local bam file -> copy it to be reachable from dragen:
			$n_bam_dragen = $dragen_input_folder.basename($n_bam);
			$parser->copyFile($n_bam, $n_bam_dragen);
		}
		
		// create cmd for vc_dragen_somatic.php
		$args = array();
		$args[] = "-t_bam ".$t_bam_dragen;
		$args[] = "-out ".$dragen_output_vcf;
		$args[] = "-out_sv ".$dragen_output_svs;
		$args[] = "-build ".$sys['build'];
		$args[] = "--log ".$dragen_log_file;
		
		//TODO again in Dragen 4.4 - Test sample: DNA2506018A1_01-DNA2504724A1_01 or DNA2505658A1_01-DNA2504802A1_01
		// $dragen_normal_vcf = $n_folder."/dragen_variant_calls/{$n_id}_dragen.vcf.gz";
		// if ($sys['type'] == "WGS" && is_file($dragen_normal_vcf))
		// {
			//calc dragen CNVs for WGS samples to compare results to clincnv
			// $args[] = "-out_cnv ".$dragen_output_cnvs;
			// $args[] = "-normal_snvs ".$dragen_normal_vcf;
		// }
		
		$args[] = "-n_bam ".$n_bam_dragen;
		
		if ($sys['type'] != "WGS" && $sys['type'] != "WGS (shallow)")
		{
			$args[] = "-is_targeted";
		}
		
		$cmd = "php ".realpath(repository_basedir())."/src/Tools/vc_dragen_somatic.php ".implode(" ", $args);
		
		// submit GridEngine job to dragen queue
		$dragen_queues = explode(",", get_path("queues_dragen"));
		$sge_args = array();
		$sge_args[] = "-V";
		$sge_args[] = "-b y"; // treat as binary
		$sge_args[] = "-wd $dragen_output_folder";
		$sge_args[] = "-m n"; // switch off messages
		$sge_args[] = "-e ".get_path("dragen_log")."/$sge_logfile.err"; // stderr
		$sge_args[] = "-o ".get_path("dragen_log")."/$sge_logfile.out"; // stdout
		$sge_args[] = "-q ".implode(",", $dragen_queues); // define queue
		$sge_args[] = "-N megSAP_DRAGEN_{$t_id}_{$n_id}"; // set name
		$qsub_command_args = implode(" ", $sge_args)." ".$cmd;

		// log sge command
		$parser->log("SGE command:\tqsub {$qsub_command_args}");

		// run qsub as user bioinf
		list($stdout, $stderr) = $parser->exec("qsub", $qsub_command_args);
		$sge_id = explode(" ", $stdout[0])[2];

		// check if submission was successful
		if ($sge_id<=0) 
		{
			trigger_error("SGE command failed:\nqsub {$qsub_command_args}\nSTDOUT:\n".implode("\n", $stdout)."\nSTDERR:\n".implode("\n", $stderr), E_USER_ERROR);
		}

		// wait for job to finish
		do 
		{
			// wait
			sleep($sge_update_interval);

			// check if job is still running
			list($stdout) = exec2("qstat -u '*' | egrep '^\s*{$sge_id}\s+' 2>&1", false);
			$finished = trim(implode("", $stdout))=="";

			// log running state
			if (!$finished)
			{
				$state = explode(" ", preg_replace('/\s+/', ' ', $stdout[0]))[4];
				trigger_error("SGE job $sge_id still queued/running (state: {$state}).", E_USER_NOTICE);
			}
		}
		while (!$finished);
		
		trigger_error("SGE job $sge_id finished.", E_USER_NOTICE);
		
		//parse SGE stdout file to determine if mapping was successful
		if (!(file_exists(get_path("dragen_log")."/$sge_logfile.out") && (get_path("dragen_log")."/$sge_logfile.err")))
		{
			trigger_error("Cannot find log files of DRAGEN mapping SGE job at '".get_path("dragen_log")."/$sge_logfile.out'!", E_USER_ERROR);
		}
		$sge_stdout = file(get_path("dragen_log")."/$sge_logfile.out", FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
		$sge_stderr = file(get_path("dragen_log")."/$sge_logfile.err", FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);

		//copy DRAGEN log to current mapping log and delete it
		if (!file_exists($dragen_log_file))
		{
			trigger_error("Cannot find DRAGEN log file '$dragen_log_file'!", E_USER_ERROR);
		}
		$parser->log("DRAGEN somatic calling log: ", file($dragen_log_file));
		unlink($dragen_log_file);

		if (end($sge_stdout)=="DRAGEN successfully finished!")
		{
			trigger_error("SGE job $sge_id successfully finished with exit status 0.", E_USER_NOTICE);
		}
		else
		{
			// write dragen log stdout and stderr to log:
			$parser->log("sge stdout:", $sge_stdout);
			$parser->log("sge stderr:", $sge_stderr);
			trigger_error("SGE job $sge_id failed!", E_USER_ERROR);
		}
		
		//copy small variant calls from Dragen
		$dragen_call_folder = $out_folder."/dragen_variant_calls/";
		if (!file_exists($dragen_call_folder))
		{
			if (!mkdir($dragen_call_folder))
			{
				trigger_error("Could not create DRAGEN variant calls folder: ".$dragen_call_folder, E_USER_ERROR);
			}
		}
		//copy out files
		$parser->moveFile($dragen_output_vcf, $dragen_call_folder.basename($dragen_output_vcf));
		$parser->moveFile($dragen_output_vcf.".tbi", $dragen_call_folder.basename($dragen_output_vcf).".tbi");
		
		if (file_exists($dragen_output_msi))
		{
			$parser->moveFile($dragen_output_msi, $dragen_call_folder.basename($dragen_output_msi));
		}
		
		if (file_exists($dragen_output_cnvs))
		{
			$parser->moveFile($dragen_output_cnvs, $dragen_call_folder.basename($dragen_output_cnvs));
			$parser->moveFile($dragen_output_cnvs.".tbi", $dragen_call_folder.basename($dragen_output_cnvs).".tbi");
		}
		
		//filter dragen vcf
		$args = array();
		$args[] = "-in ".$dragen_call_folder.basename($dragen_output_vcf);
		$args[] = "-tumor_name {$t_id}";
		$args[] = "-normal_name {$n_id}";
		$args[] = "-out {$variants}";
		$args[] = "-build ".$sys['build'];
		$args[] = "-target ".$roi;
		$parser->execTool("Tools/an_filter_dragen_somatic.php", implode(" ", $args));
		
		
		if (is_file($dragen_output_svs))
		{
			//copy svs
			$parser->moveFile($dragen_output_svs, $dragen_call_folder.basename($dragen_output_svs));
			$parser->moveFile($dragen_output_svs.".tbi", $dragen_call_folder.basename($dragen_output_svs).".tbi");
				
			$parser->copyFile($dragen_call_folder.basename($dragen_output_svs), $manta_sv);
			$parser->copyFile($dragen_call_folder.basename($dragen_output_svs).".tbi", $manta_sv.".tbi");
			
			$parser->execApptainer("ngs-bits", "VcfToBedpe", "-in $manta_sv -out $manta_sv_bedpe", [$manta_sv], [dirname($manta_sv_bedpe)]);
			$parser->execTool("Tools/bedpe2somatic.php", "-in $manta_sv_bedpe -out $manta_sv_bedpe -tid $t_id -nid $n_id");
			
			if( db_is_enabled("NGSD") )
			{
				$parser->execApptainer("ngs-bits", "BedpeGeneAnnotation", "-in $manta_sv_bedpe -out $manta_sv_bedpe -add_simple_gene_names", [$manta_sv_bedpe]);
			}
		}
		
		#remove copied files if exist:
		if (contains($t_bam_dragen, $dragen_input_folder))
		{
			#if file was copied into the dragen input folder delete it now.
			unlink($t_bam_dragen);
		}
		if (contains($n_bam_dragen, $dragen_input_folder))
		{
			#if file was copied into the dragen input folder delete it now.
			unlink($n_bam_dragen);
		}
		
	}
	elseif($use_deepsomatic)
	{
		$args = [];

		if ($sys['type'] === "WGS")	$args[] = "-model_type WGS";
		else $args[] = "-model_type WES";

		$args[] = "-bam_tumor ".$t_bam;
		$args[] = "-bam_normal ".$n_bam;
		$args[] = "-out ".$variants;
		$args[] = "-build ".$sys['build'];
		$args[] = "-threads ".$threads;
		$args[] = "-tumor_id {$t_id}";
		$args[] = "-normal_id {$n_id}";
		$args[] = "-default";

		if (!empty($roi))
		{
			$args[] = "-target {$roi}";
		}
		#$args[] = "-allow_empty_examples";

		$parser->execTool("Tools/vc_deepsomatic.php", implode(" ", $args));
	}
	else //Strelka calling
	{
		$args_strelka = [
			"-t_bam {$t_bam}",
			"-n_bam {$n_bam}",
			"-out {$variants}",
			"-build ".$sys['build'],
			"-threads {$threads}"
		];
		if (!empty($roi))
		{
			$args_strelka[] = "-target {$roi}";
		}
		if ($sys['type'] === "WGS")
		{
			$args_strelka[] = "-wgs";
		}
		if (is_file($manta_indels))
		{
			$args_strelka[] = "-smallIndels {$manta_indels}";
		}
		$parser->execTool("Tools/vc_strelka2.php", implode(" ", $args_strelka));
	}

	//add somatic BAF file
	$variants_germline_vcf = "{$n_basename}_var.vcf.gz";
	if (file_exists($variants_germline_vcf))
	{
		$baf_args = [
			"-bam_t {$t_bam}",
			"-bam_n {$n_bam}",
			"-vcf {$variants_germline_vcf}",
			"-out {$ballele}",
			"-build ".$sys['build']
		];
		
		if ($sys['type'] === "WGS")
		{
			$baf_args[] = "-downsample 100";
		}
		
		$parser->execTool("Tools/baf_somatic.php", implode(" ", $baf_args));
	}
	else
	{
		trigger_error("Cannot create BAF file as normal VCF missing: {$variants_germline_vcf}", E_USER_NOTICE);
	}
	
	if (file_exists($variants) && !$skip_signatures)
	{
		$snv_signatures_out = $out_folder."/snv_signatures/";
		$tmp_variants = $parser->tempFile(".vcf", "snv_signatures_");
		$parser->execApptainer("htslib", "bgzip", "-c -d -@ {$threads} $variants > {$tmp_variants}", [$variants]);
		$parser->exec("php ".repository_basedir()."/src/Tools/extract_signatures.php", "-in {$tmp_variants} -mode snv -out {$snv_signatures_out} -reference GRCh38 -threads {$threads}", true);
	}
}

//detection of viral DNA
if (in_array("vi", $steps))
{
	$t_bam_dedup = "{$t_basename}_before_dedup.bam";
	$dedup_used = file_exists($t_bam_dedup);
	$vc_viral_args = [
		"-in ".($dedup_used ? $t_bam_dedup : $t_bam),
		"-viral_bam {$t_basename}_viral.bam",
		"-viral_bam_raw {$t_basename}_viral_before_dedup.bam",
		"-viral_cov {$t_basename}_viral.tsv",
		"-viral_chrs chrNC_007605",
		"-build_viral somatic_viral",
		"-avg_target_cov ".get_qcml_value("{$t_basename}_stats_map.qcML", "QC:2000025"),
		"-threads ".$threads,
	];
	if ($dedup_used) $vc_viral_args[] = "-barcode_correction";
	$parser->execTool("Tools/vc_viral_load.php", implode(" ", $vc_viral_args));
}

//CNV calling
$som_clincnv = $full_prefix . "_clincnv.tsv"; //ClinCNV output file
if(in_array("cn",$steps))
{
	//for differential analysis it is necessary that both tumor and normal sample have the same processing system.
	if ($sys["name_short"] != $n_sys["name_short"])
	{
		trigger_error("Tumor and normal sample were sequenced with different processing systems - CNVs cannot be calculated!",E_USER_ERROR);
	}
	
	// copy number variant calling
	$tmp_folder = $parser->tempFolder();
	//directory with reference coverage files and directory with reference off-target coverage files
	$ref_folder_t = get_path("data_folder")."/coverage/".$sys['name_short']."-tumor";
	$ref_folder_t_off_target = $ref_folder_t . "_off_target";
	
	$ref_folder_n = get_path("data_folder")."/coverage/".$sys['name_short']; //same as $n_sys['name_short"] if available
	$ref_folder_n_off_target = $ref_folder_n . "_off_target";
	
	create_directory($ref_folder_t);
	create_directory($ref_folder_t_off_target);
	create_directory($ref_folder_n);
	create_directory($ref_folder_n_off_target);
	
	//Generate file with off-target region
	$off_target_bed = get_path("data_folder")."/coverage/off_target_beds/".$sys['name_short'].".bed";
	if(!file_exists($off_target_bed)) create_off_target_bed_file($off_target_bed,$sys['target_file'], $ref_genome);
	
	$target_bed = "";
	
	//get correct coverage files, target bed and off-target bed based on system type:
	if ($sys["type"] != "WGS" && $sys['type'] != "WGS (shallow)")
	{
		$target_bed = $ref_folder_t."/roi_annotated.bed";
		if (!file_exists($target_bed))
		{
			$pipeline = [
					["", $parser->execApptainer("ngs-bits", "BedAnnotateGC", "-in ".realpath($sys['target_file'])." -clear -ref {$ref_genome}", [$sys['target_file'], $ref_genome], [], true)],
					["", $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-out {$target_bed}", [], [dirname($target_bed)], true)],
				];
			$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
		}
	}
	else // WGS sample
	{
		$bin_size = get_path("cnv_bin_size_wgs");
		
		//target bed
		$target_bed = $ref_folder_t."/bins{$bin_size}.bed"; // is outside of bins folder
		if (!file_exists($target_bed))
		{
			$pipeline = [
					["", $parser->execApptainer("ngs-bits", "BedChunk", "-in ".realpath($sys['target_file'])." -n {$bin_size}", [$sys['target_file']], [], true)],
					["", $parser->execApptainer("ngs-bits", "BedAnnotateGC", "-clear -ref {$ref_genome}", [$ref_genome], [], true)],
					["", $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-out {$target_bed}", [], [dirname($target_bed)], true)]
				];
			$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
		}
		
		//coverage files:
		
		$bin_folder_t = "{$ref_folder_t}/bins{$bin_size}/";
		$bin_folder_n = "{$ref_folder_n}/bins{$bin_size}/";
		create_directory($bin_folder_t);
		create_directory($bin_folder_n);
		$ref_folder_t = $bin_folder_t;
		$ref_folder_n = $bin_folder_n;
	}
	
	
	/***************************************************
	 * GENERATE AND COPY COVERAGE FILES TO DATA FOLDER *
	 ***************************************************/
	// coverage for tumor sample + off-target files
	$t_cov = "{$tmp_folder}/{$t_id}.cov";
	$t_cov_off_target = "{$tmp_folder}/{$t_id}_off_target.cov";
	$ref_file_t = "{$ref_folder_t}/{$t_id}.cov";
	$ref_file_t_off_target = "{$ref_folder_t_off_target}/{$t_id}.cov";
	
	$parser->execApptainer("ngs-bits", "BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam $t_bam -in $target_bed -out $t_cov -threads {$threads} -ref {$ref_genome}", [$t_bam, $target_bed, $ref_genome]);
	$parser->execApptainer("ngs-bits", "BedCoverage", "-clear -min_mapq 10 -decimals 4 -in $off_target_bed -bam $t_bam -out $t_cov_off_target -threads {$threads} -ref {$ref_genome}", [$off_target_bed, $t_bam, $ref_genome]);

	
	//copy tumor sample coverage file to reference folder (has to be done before ClinCNV call to avoid analyzing the same sample twice)
	if (db_is_enabled("NGSD") && is_valid_ref_tumor_sample_for_cnv_analysis($t_id))
	{
		$parser->log("Moving tumor sample coverage file to reference folder...");
		$parser->copyFile($t_cov, $ref_file_t); 
		$t_cov = $ref_file_t;

		$parser->copyFile($t_cov_off_target,$ref_file_t_off_target);
		$t_cov_off_target = $ref_file_t_off_target;
	}
	
	// coverage file for normal sample
	$n_cov = "{$ref_folder_n}/{$n_id}.cov.gz";
	if (!file_exists($n_cov))
	{
		$cov_tmp_unzipped = $tmp_folder."/{$n_id}.cov";
		$parser->execApptainer("ngs-bits", "BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam $n_bam -in $target_bed -out $cov_tmp_unzipped -threads {$threads} -ref {$ref_genome}", [$n_bam, $target_bed, $ref_genome]);
		$parser->exec("gzip", "-9 {$cov_tmp_unzipped}");
		
		if (db_is_enabled("NGSD") && is_valid_ref_sample_for_cnv_analysis($n_id))
		{
			$parser->log("Moving normal sample coverage file to reference folder...");
			$parser->copyFile($cov_tmp_unzipped.".gz", $n_cov);
		}
		else
		{
			$n_cov = $cov_tmp_unzipped.".gz";
		}
	}
	
	// coverage file for normal sample (off-target)
	$n_cov_off_target = "{$tmp_folder}/{$n_id}_off_target.cov";
	$ref_file_n_off_target = "{$ref_folder_n_off_target}/{$n_id}.cov";
	$parser->execApptainer("ngs-bits", "BedCoverage", "-clear -min_mapq 10 -decimals 4 -in $off_target_bed -bam $n_bam -out $n_cov_off_target -threads {$threads} -ref {$ref_genome}", [$off_target_bed, $n_bam, $ref_genome]);
	if (db_is_enabled("NGSD") && is_valid_ref_sample_for_cnv_analysis($n_id))
	{
		$parser->log("Moving normal sample off-target coverage file to reference folder...");
		$parser->copyFile($n_cov_off_target, $ref_file_n_off_target);
		$n_cov_off_target = $ref_file_n_off_target;
	}

	//append tumor-normal IDs to list with tumor normal IDs (stored in same folder as tumor coverage files)
	$t_n_list_file = $ref_folder_t . "/" . "list_tid-nid.csv";

	if (!file_exists($t_n_list_file))
	{
		$header = "##THIS FILE CONTAINS TUMOR AND NORMAL IDS OF PROCESSING SYSTEM ".$sys['name_short']."\n";
		$header .= "#tumor_id,normal_id\n";
		file_put_contents($t_n_list_file,$header);
	}
	
	//use temporary list file if n or t cov files are not valid
	if(!db_is_enabled("NGSD") || !is_valid_ref_sample_for_cnv_analysis($n_id) || !is_valid_ref_tumor_sample_for_cnv_analysis($t_id))
	{
		$tmp_file_name = $parser->tempFile(".csv");
		$parser->copyFile($t_n_list_file,$tmp_file_name);
		$t_n_list_file = $tmp_file_name;
	}
	
	$baf_folder = get_path("data_folder")."/coverage/". $sys['name_short']."_bafs";
	create_directory($baf_folder);
	if(db_is_enabled("NGSD"))
	{
		$db_conn = DB::getInstance("NGSD");
		//Create BAF file for each sample with the same processing system if not existing
		complement_baf_folder($t_id, $n_id, $t_n_list_file,$baf_folder,$db_conn, $sys['build']);
	}
	else
	{
		$normal_gsvar = "{$n_basename}.GSvar";
		$parser->execTool("Auxilary/create_baf_file.php", "-gsvar $normal_gsvar -bam $n_bam -build ".$sys['build']." -out_file {$baf_folder}/{$n_id}.tsv");
	}

	//Skip CNV Calling if there are less than specified tumor-normal coverage pairs
	if(count(file($t_n_list_file)) > $min_cov_files)
	{
		/*******************
		 * EXECUTE CLINCNV *
		 *******************/
		$cohort_folder = get_path("clincnv_cohorts")."/". $sys['name_short'];
		if(!file_exists($cohort_folder)) mkdir($cohort_folder, 0777);
			
		$args_clincnv = [
		"-t_id", $t_id,
		"-n_id", $n_id,
		"-t_cov", $t_cov,
		"-n_cov", $n_cov,
		"-out", $som_clincnv,
		"-cov_pairs",$t_n_list_file,
		"-system", $system,
		"-bed", $target_bed,
		"-baf_folder", $baf_folder,
		"-cohort_folder", $cohort_folder,
		];
		
		if ($sys["type"] != "WGS" && $sys["type"] != "WGS (shallow)")
		{
			$args_clincnv[] = "-bed_off {$off_target_bed}";
			$args_clincnv[] = "-t_cov_off {$t_cov_off_target}";
			$args_clincnv[] = "-n_cov_off {$n_cov_off_target}";
			$args_clincnv[] = "-cov_folder_t_off {$ref_folder_t_off_target}";
			$args_clincnv[] = "-cov_folder_n_off {$ref_folder_n_off_target}";
			$args_clincnv[] = "-threads {$threads}";
		}

		if($sys['type'] == "WES")
		{
			$args_clincnv[] = "-lengthS 9";
			$args_clincnv[] = "-scoreS 200";
			$args_clincnv[] = "-filterStep 2";
		}
		
		if($sys['type'] == "WGS")
		{
			$args_clincnv[] = "-lengthS 9";
			$args_clincnv[] = "-scoreS 1500";
			$args_clincnv[] = "-filterStep 2";
			$args_clincnv[] = "-clonePenalty 10000";
			//for WGS samples use at max 2 threads as it uses a lot of RAM
			$args_clincnv[] = "-threads ".min(2, $threads);
		}
		
		if(isset($cnv_baseline_pos))
		{
			$args_clincnv[] = "-guide_baseline $cnv_baseline_pos";
		}
		
		if(isset($cnv_wgs_purity_step))
		{
			$args_clincnv[] = "-purityStep $cnv_wgs_purity_step";
		}
		
		$parser->execTool("Tools/vc_clincnv_somatic.php",implode(" ",$args_clincnv));
		
		//Annotate cytoband and data from network of cancer genes
		$parser->execTool("Tools/an_somatic_cnvs.php","-cnv_in $som_clincnv -out $som_clincnv -include_ncg -include_cytoband");
	}
	else
	{
		trigger_error("Not enough reference tumor-normal coverage files for processing system {$system} found. Skipping CNV calling.\n", E_USER_NOTICE);
	}
	
	//calculate HRD based on clincnv.tsv file
	if (file_exists($som_clincnv))
	{
		if (!$skip_HRD)
		{
			$parser->execTool("Tools/an_scarHRD.php" , "-cnvs {$som_clincnv} -tumor {$t_id} -normal {$n_id} -out_folder {$out_folder}");
		}
		
		if(!$skip_signatures)
		{
			$cnv_signatures_out = $out_folder."/cnv_signatures/";
			$parser->exec("php ".repository_basedir()."/src/Tools/extract_signatures.php", "-in {$som_clincnv} -mode cnv -out {$cnv_signatures_out} -reference GRCh38 -threads {$threads}", true);
		}
	}
}

//annotation
$variants_annotated = $full_prefix . "_var_annotated.vcf.gz";	// annotated variants
$variants_gsvar     = $full_prefix . ".GSvar";					// GSvar variants
$somaticqc          = $full_prefix . "_stats_som.qcML";			// SomaticQC qcML
$cfdna_folder       = $full_prefix . "_cfDNA_candidates";       // folder containing cfDNA monitoring variants 
if (in_array("an", $steps))
{
	check_genome_build($variants, $sys['build']);
	
	// annotate vcf (in temp folder)
	$tmp_folder1 = $parser->tempFolder();
	$tmp_vcf = "{$tmp_folder1}/{$prefix}_var_annotated.vcf.gz";
	$parser->execTool("Tools/annotate.php", "-out_name $prefix -out_folder $tmp_folder1 -system $system -vcf $variants -somatic -threads $threads");

	// run somatic QC
	$links = array_filter([
		"{$t_basename}_stats_fastq.qcML",
		"{$t_basename}_stats_map.qcML",
		"{$n_basename}_stats_fastq.qcML",
		"{$n_basename}_stats_map.qcML"
	], "file_exists");
	
	$in_files = [
		$t_bam,
		$n_bam,
		$roi,
		repository_basedir()."/data/gene_lists",
		$ref_genome
	];

	$args_somaticqc = [
		"-tumor_bam", $t_bam,
		"-normal_bam", $n_bam,
		"-somatic_vcf", $tmp_vcf,
		"-target_bed", realpath($roi),
		"-target_exons", repository_basedir()."/data/gene_lists/gene_exons.bed", //file containing all human exons to determine exonic variants in TMB calculation
		"-blacklist", repository_basedir() ."/data/gene_lists/somatic_tmb_blacklist.bed", //Blacklisted genes that are not included in TMB calculation (e.g. HLA-A and HLA-B)
		"-tsg_bed", repository_basedir() ."/data/gene_lists/somatic_tmb_tsg.bed", //TSG genes whose mutations are treated specially in TMB calculation
		"-ref", $ref_genome,
		"-out", $somaticqc,
		"-build", ngsbits_build($sys['build'])
	];
	if (!empty($links))
	{
		$args_somaticqc[] = "-links";
		$args_somaticqc[] = implode(" ", $links);
		$in_files = array_merge($in_files, $links);
	}
	$parser->execApptainer("ngs-bits", "SomaticQC", implode(" ", $args_somaticqc), $in_files, [dirname($somaticqc)]);

	//add sample info to VCF header
	$s = Matrix::fromTSV($tmp_vcf);
	$comments = $s->getComments();
	$comments[] = gsvar_sample_header($t_id, array("IsTumor" => "yes"), "#", "");
	$comments[] = gsvar_sample_header($n_id, array("IsTumor" => "no"), "#", "");
	$s->setComments(sort_vcf_comments($comments));
	$s->toTSV($tmp_vcf);

	// zip and index vcf file
	$parser->execApptainer("htslib", "bgzip", "-c $tmp_vcf > $variants_annotated", [], [dirname($variants_annotated)]);
	$parser->execApptainer("htslib", "tabix", "-f -p vcf $variants_annotated", [], [dirname($variants_annotated)]);

	// convert vcf to GSvar
	$args = array("-in $tmp_vcf", "-out $variants_gsvar", "-t_col $t_id", "-n_col $n_id");
	$parser->execTool("Tools/vcf2gsvar_somatic.php", implode(" ", $args));
	
	//Annotate data from network of cancer genes
	$parser->execTool("Tools/an_somatic_gsvar.php" , "-gsvar_in $variants_gsvar -out $variants_gsvar -include_ncg");

	//Determine cfDNA monitoring candidates
	//remove previous calls
	if (file_exists($cfdna_folder)) exec2("rm -r $cfdna_folder"); 

	// call variant selection in umiVar container
	$params = array();
	$params[] = "-v $tmp_vcf";
	$params[] = "-g $variants_gsvar";
	$params[] = "-o $cfdna_folder";
	$params[] = "-r $ref_genome";
	$params[] = "-n 65"; //number of variants to select
	$params[] = "-i"; // ignore INDELS
	$params[] = "-s $n_id";
	$parser->execApptainer("umiVar", "select_monitoring_variants.py", implode(" ", $params), [$variants_gsvar, $ref_genome], [$out_folder]);
}

//MSI calling
$msi_o_file = $full_prefix . "_msi.tsv";						//MSI
if (in_array("msi", $steps))
{
	//temp folder so other output (_dis, _germline and _somatic files) is written to tmp and automatically deleted even if the tool crashes
	$msi_tmp_folder = $parser->tempFolder("sp_detect_msi_");
	$msi_tmp_file = $msi_tmp_folder."/msi_tmp_out.tsv";
	$msi_folder = get_path("data_folder")."/dbs/msisensor-pro/";
	
	if (!file_exists($msi_folder)) $parser->exec("mkdir", "-p {$msi_folder}");
	
	$msi_ref = $msi_folder."/msisensor_references_".$n_sys['build'].get_path("container_msisensor-pro").".site";
	$parameters = "-n_bam $n_bam -t_bam $t_bam -msi_ref $msi_ref -ref $ref_genome -threads $threads -out " .$msi_tmp_file. " -build ".$n_sys['build'];
	
	$parser->execTool("Tools/detect_msi.php", $parameters);
	$parser->copyFile($msi_tmp_file, $msi_o_file);
}

if (in_array("an_rna", $steps))
{
	$args = array();
	$args[] = "-t_bam $t_bam";
	$args[] = "-full_prefix $full_prefix";
	$args[] = "-steps ".count($steps);
	if (file_exists($system)) $args[] = "-system $system";
	if (file_exists($t_rna_bam)) $args[] = "-t_rna_bam $t_rna_bam";
	if ($rna_ref_tissue != "") $args[] = "-rna_ref_tissue $rna_ref_tissue";
	if ($skip_correlation) $args[] = "-skip_correlation";

	$parser->execTool("Tools/an_somatic_rna.php", implode(" ", $args));
}

//Collect QC terms if necessary
$qc_other = $full_prefix."_stats_other.qcML";
if (in_array("vc", $steps) || in_array("vi", $steps) || in_array("msi", $steps) || in_array("cn", $steps) || in_array("db", $steps))
{
	$parser->execTool("Tools/create_qcml.php", "-full_prefix {$full_prefix} -viral_tsv {$t_basename}_viral.tsv");
}

//DB import
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	
	$t_info = get_processed_sample_info($db_conn, $t_id, false);
	$n_info = get_processed_sample_info($db_conn, $n_id, false);
		
	if (is_null($t_info) || (is_null($n_info)))
	{
		trigger_error("No database import since no valid processing ID (T: {$t_id} /N: {$n_id})", E_USER_WARNING);
	}
	else
	{
		// import qcML files
		$qcmls = implode(" ", array_filter([
			"{$t_basename}_stats_fastq.qcML",
			"{$t_basename}_stats_map.qcML",
			$somaticqc,
			$qc_other
		], "file_exists"));
		$parser->execApptainer("ngs-bits", "NGSDImportSampleQC", "-ps $t_id -files $qcmls -force", ["{$t_basename}_stats_fastq.qcML", "{$t_basename}_stats_map.qcML", $somaticqc, $qc_other]);

		// check tumor/normal flag
		if (!$t_info['is_tumor'])
		{
			trigger_error("Tumor sample $t_id is not flagged as tumor in NGSD!", E_USER_WARNING);
		}

		// check sex using control sample
		$parser->execTool("Tools/db_check_gender.php", "-in $n_bam -pid $n_id -sry_cov 30");

		// check tumor/normal flag
		if ($n_info['is_tumor'])
		{
			trigger_error("Normal sample $n_id is flagged as tumor in NGSD!", E_USER_WARNING);
		}

		// update normal sample entry for tumor
		if (updateNormalSample($db_conn, $t_id, $n_id))
		{
			trigger_error("Updated normal sample ($n_id) for tumor ($t_id) in NGSD.", E_USER_NOTICE);
		}
		
		// import sample relation
		$s_id_t = $t_info['s_id'];
		$s_id_n = $n_info['s_id'];
		$db_conn->executeStmt("INSERT IGNORE INTO `sample_relations`(`sample1_id`, `relation`, `sample2_id`) VALUES ({$s_id_t},'tumor-normal',{$s_id_n})");
		
		// import variants into NGSD
		$args = ["-t_ps {$t_id}", "-n_ps {$n_id}", "-force"];
		$binds = [];
		if (file_exists($variants_gsvar))
		{
			check_genome_build($variants_gsvar, $sys['build']);
			$args[] = "-var {$variants_gsvar}";
			$binds[] = $variants_gsvar;
		}
		if(file_exists($som_clincnv))
		{
			$args[] = "-cnv {$som_clincnv}";
			$binds[] = $som_clincnv;
		}			
		if(file_exists($manta_sv_bedpe))
		{
			$args[] = "-sv {$manta_sv_bedpe}";
			$binds[] = $manta_sv_bedpe;
		}
		if (count($binds)>0)
		{
			$parser->execApptainer("ngs-bits", "NGSDAddVariantsSomatic", implode(" ", $args), $binds);
		}	
		
		//add secondary analysis (if missing)
		$parser->execTool("Tools/db_import_secondary_analysis.php", "-type 'somatic' -gsvar {$variants_gsvar}");
	}
}
?>
