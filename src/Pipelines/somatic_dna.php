<?php

/**
 * @page somatic_dna
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_dna", "Differential analysis of tumor and normal DNA samples.");

//mandatory
$parser->addInfile("t_bam", "Tumor sample BAM file.", false);
$parser->addString("out_folder", "Output folder.", false);

//optional
$parser->addInfile("n_bam", "Normal sample BAM file.", true);
$parser->addInfile("t_rna_bam", "Tumor RNA sample BAM file.", true);
$parser->addString("prefix", "Output file prefix.", true, "somatic");

$steps_all = array("vc", "vi", "cn", "an", "msi", "an_rna", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\n" .
	"vc=variant calling, an=annotation,\n" .
	"cn=copy-number analysis, msi=microsatellite analysis,\n".
	"an_rna=annotate data from somatic RNA files,\n".
	"vi=virus detection, db=database import",
	true, "vc,cn,an,msi,vi,an_rna,db");

$parser->addInfile("system",  "Processing system file used for tumor DNA sample (resolved from NGSD via tumor BAM by default).", true);
$parser->addInfile("n_system",  "Processing system file used for normal DNA sample (resolved from NGSD via normal BAM by default).", true);
$parser->addFlag("skip_contamination_check", "Skips check of female tumor sample for male SRY DNA.");
$parser->addFlag("skip_correlation", "Skip sample correlation check.");
$parser->addFlag("skip_low_cov", "Skip low coverage statistics.");
$parser->addFlag("skip_signatures", "Skip calculation of mutational signatures.");
$parser->addFlag("skip_HRD", "Skip calculation HRD.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
$parser->addFlag("use_dragen", "Use Illumina dragen for somatic variant calling.");
//default cut-offs
$parser->addFloat("min_af", "Allele frequency detection limit (for tumor-only calling only).", true, 0.05);
$parser->addFloat("min_correlation", "Minimum correlation for tumor/normal pair.", true, 0.8);
$parser->addFloat("min_depth_t", "Tumor sample coverage cut-off for low coverage statistics.", true, 60);
$parser->addFloat("min_depth_n", "Normal sample coverage cut-off for low coverage statistics.", true, 60);
$parser->addInt("min_cov_files", "Minimum number of required tumor-normal pairs for CNV calling.", true, 7);
$parser->addString("cnv_baseline_pos","baseline region for ClinCNV, format e.g. chr1:12-12532",true);
$parser->addString("rna_ref_tissue", "Reference data for RNA annotation", true);
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));


###################################### AUXILARY FUNCTIONS ######################################

//Creates BAF file from "$gsvar" and "bam" file and writes content to "out_file"
function create_baf_file($gsvar, $bam, $out_file, $ref_genome, &$error = false)
{
	global $parser;

	if(!file_exists($gsvar) || !file_exists($bam)) 
	{
		trigger_error("Could not create BAF file {$out_file}, no GSvar or BAM file available.", E_USER_WARNING);
		$error = true;
		return;
	}
	//Abort if out_file exists to prevent interference with other jobs
	if(file_exists($out_file)) return;

	$tmp_out = $parser->tempFile(".tsv");
	exec2(get_path("ngs-bits")."/VariantAnnotateFrequency -in {$gsvar} -bam {$bam} -depth -out {$tmp_out} -ref {$ref_genome}", true);

	$in_handle  = fopen2($tmp_out,"r");
	$out_handle = fopen2($out_file,"w");

	while(!feof($in_handle))
	{

		$line = trim(fgets($in_handle));
		if(starts_with($line,"#")) continue;
		if(empty($line)) continue;
		
		$parts = explode("\t",$line);
		list($chr,$start,$end,$ref,$alt) = $parts;
		if(strlen($ref) != 1 || strlen($alt) != 1 ) continue; //Skip INDELs

		$freq = $parts[count($parts)-2]; //freq annotation ist 2nd last col
		$depth= $parts[count($parts)-1]; //depth annotation is last col
		fputs($out_handle,"{$chr}\t{$start}\t{$end}\t{$chr}_{$start}\t{$freq}\t{$depth}\n");
	}
	fclose($in_handle);
	fclose($out_handle);
}

//Checks $baf_folder for missing B-AF files and creates them if neccessary
function complement_baf_folder($t_n_id_file,$baf_folder,&$db_conn, $ref_genome)
{
	$ids = file($t_n_id_file);
	foreach($ids as $line)
	{
		if(starts_with($line,"#")) continue;
		list($tid,$nid) = explode(",",trim($line));
		if(!file_exists("{$baf_folder}/{$nid}.tsv"))
		{
			$ninfo = get_processed_sample_info($db_conn,$nid);
			$n_gsvar = $ninfo["ps_folder"] ."/{$nid}.GSvar";
			$n_bam = $ninfo["ps_bam"];
			create_baf_file($n_gsvar,$n_bam,"{$baf_folder}/{$nid}.tsv", $ref_genome);
		}
		if(!file_exists("{$baf_folder}/{$tid}.tsv"))
		{
			$ninfo = get_processed_sample_info($db_conn,$nid);
			$n_gsvar = $ninfo["ps_folder"] ."/{$nid}.GSvar";
			$tinfo = get_processed_sample_info($db_conn,$tid);
			$t_bam = $tinfo["ps_bam"];
			create_baf_file($n_gsvar,$t_bam,"{$baf_folder}/{$tid}.tsv", $ref_genome);
		}
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
	
	if (!isset($n_bam))
	{
		trigger_error("Dragon analysis currently not supported for single sample analysis!", E_USER_ERROR);
	}
}


###################################### SCRIPT START ######################################
if (!file_exists($out_folder))
{
	exec2("mkdir -p $out_folder");
}
if ($parser->getLogFile() == "") $parser->setLogFile($out_folder."/somatic_dna_".date("YmdHis").".log");

//output prefix
$full_prefix = "{$out_folder}/{$prefix}";

//IDs, system and target region
$t_id = basename2($t_bam);
$t_basename = dirname($t_bam)."/".$t_id;
$sys = load_system($system, $t_id);
$roi = $sys["target_file"];
$ref_genome = genome_fasta($sys['build']);

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);
}

//make sure it is a BAM (MSIsensor does not work on CRAM)
check_genome_build($t_bam, $sys['build']);
$t_bam = convert_to_bam_if_cram($t_bam, $parser, $sys['build'], $threads);

//normal sample data (if not single sample analysis)
$single_sample = !isset($n_bam);
if (!$single_sample)
{
	$n_id = basename2($n_bam);
	$n_basename = dirname($n_bam)."/".$n_id;
	$n_sys = load_system($n_system, $n_id);
	$ref_folder_n = get_path("data_folder")."/coverage/".$n_sys['name_short'];
	
	//make sure it is a BAM (MSIsensor does not work on CRAM)
	check_genome_build($n_bam, $n_sys['build']);
	$n_bam = convert_to_bam_if_cram($n_bam, $parser, $sys['build'], $threads);
	
	//Check whether both samples have same processing system
	if($roi != $n_sys["target_file"])
	{
		#test that tumor target is a subset of normal target
		exec(get_path("ngs-bits")."BedSubtract -in ".$roi." -in2 ".$n_sys["target_file"], $output, $return_var);
		
		foreach ($output as $line)
		{
			if ($line == "" || starts_with($line, "#")) continue;
			
			trigger_error("Tumor sample $t_id  target region is different from, and not a subset of, the normal sample $n_id target region.",E_USER_ERROR);
		}
		
		trigger_error("Tumor sample $t_id and normal sample $n_id have different target regions.",E_USER_WARNING);
		
		
		if (in_array("cn", $steps))
		{
			trigger_error("CNVs cannot be calculated with two different target regions. Removing 'cn' step!",E_USER_WARNING);
			
			$key = array_search("cn", $steps);
			unset($steps[$key]);
		}
	}
}

//Abort if calling is requested and somatic report config exists in NGSD
if (!$single_sample && db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
	list($config_id, $config_vars_exist, $config_cnvs_exist) = somatic_report_config($db, $t_id, $n_id);
	if (in_array("vc", $steps) && $config_vars_exist)
	{
		trigger_error("Somatic report configuration with SNVs exists in NGSD! Delete somatic report configuration for reanalysis of step 'vc'.", E_USER_ERROR);
	}
	if (in_array("cn", $steps) && $config_cnvs_exist)
	{
		trigger_error("Somatic report configuration with CNVs exists in NGSD! Delete somatic report configuration for reanalysis of step 'cn'.", E_USER_ERROR);
	}
}

//sample similarity check
$bams = array_filter([$t_bam, $n_bam]);
if (count($bams) > 1)
{
    if ($skip_correlation)
    {
        trigger_error("Genotype correlation check has been disabled!", E_USER_WARNING);
    }
    else
    {
    	$args_similarity = [
    		"-in ".implode(" ", $bams),
			"-mode bam",
			"-build ".ngsbits_build($sys['build'])
		];
		if (!empty($roi))
		{
			$args_similarity[] = "-roi {$roi}";
		}
        $output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", implode(" ", $args_similarity), true);

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
}

//check SRY coverage of tumor for female samples, can be a hint for contamination with male DNA
if( db_is_enabled("NGSD") )
{
	$db = DB::getInstance("NGSD");
	$tinfo = get_processed_sample_info($db, $t_id, false);
	
	if(!is_null($tinfo) && $tinfo["gender"] == "female")
	{
		if($skip_contamination_check)
		{
			trigger_error("Skipping check of female tumor sample $t_bam for contamination with male genomic DNA.", E_USER_WARNING);
		}
		else
		{
			$out = $parser->exec(get_path("ngs-bits") . "/SampleGender",  "-in $t_bam -build ".ngsbits_build($sys['build'])." -method sry", true);
			list(,,$cov_sry) = explode("\t", $out[0][1]);

			if(is_numeric($cov_sry) && (float)$cov_sry >= 30)
			{
				trigger_error("Detected contamination of female tumor sample {$t_id} with male genomic DNA on SRY. SRY coverage is at {$cov_sry}x.", E_USER_ERROR);
			}
		}
	}
}

// Check samples are flagged correctly in NGSD
if( db_is_enabled("NGSD") && count($bams) > 1 )
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
	$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $roi -bam $t_bam -out $low_cov -cutoff $min_depth_t -threads {$threads} -ref {$ref_genome}", true);
	//combined tumor and normal low coverage files
	//normal coverage is calculated only for tumor target region
	if(!$single_sample)
	{
		$low_cov_n = $parser->tempFile("_nlowcov.bed");
		$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $roi -bam $n_bam -out $low_cov_n -cutoff $min_depth_n -threads {$threads} -ref {$ref_genome}", true);
		$parser->execPipeline([
			[get_path("ngs-bits")."BedAdd", "-in $low_cov $low_cov_n"],
			[get_path("ngs-bits")."BedMerge", "-out $low_cov"]
		], "merge low coverage BED files");
	}
	// annotate with gene names
	if (db_is_enabled("NGSD"))
	{
		$parser->exec(get_path("ngs-bits") . "BedAnnotateGenes", "-in $low_cov -extend 25 -out $low_cov", true);
	}
}

//variant calling
$manta_indels  = $full_prefix . "_manta_var_smallIndels.vcf.gz";	// small indels from manta
$manta_sv      = $full_prefix . "_manta_var_structural.vcf.gz";		// structural variants (vcf)
$manta_sv_bedpe= $full_prefix . "_manta_var_structural.bedpe"; 		// structural variants (bedpe)
$variants      = $full_prefix . "_var.vcf.gz";					// variants
$ballele       = $full_prefix . "_bafs.igv";					// B-allele frequencies
$hla_file_tumor = "{$t_basename}_hla_genotyper.tsv";
if (!$single_sample)
{
	$hla_file_normal = "{$n_basename}_hla_genotyper.tsv"; 
}
if (in_array("vc", $steps))
{
	//get HLA types
	$parser->execTool("NGS/hla_genotyper.php", "-bam $t_bam -name $t_id -out ".$hla_file_tumor);
	if(!$single_sample)
	{
		$parser->execTool("NGS/hla_genotyper.php", "-bam $n_bam -name $n_id -out " . $hla_file_normal);
	}
	
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
		if (!$single_sample)
		{
			$args_manta[] = "-bam $n_bam";
		}
		if ($sys['type'] !== "WGS") //use exome flag for non targeted / exome samples (i.e. non WGS samples)
		{
			$args_manta[] = "-exome";
		}
		if (!empty($roi))
		{
			$args_manta[] = "-target {$roi}";
		}
		$parser->execTool("NGS/vc_manta.php", implode(" ", $args_manta));
		
		exec2(get_path("ngs-bits") . "VcfToBedpe -in $manta_sv -out $manta_sv_bedpe");
		if(!$single_sample)
		{
			$parser->execTool("Tools/bedpe2somatic.php", "-in $manta_sv_bedpe -out $manta_sv_bedpe -tid $t_id -nid $n_id");
		}
		
		if( db_is_enabled("NGSD") )
		{
			$parser->exec(get_path("ngs-bits") . "BedpeGeneAnnotation", "-in $manta_sv_bedpe -out $manta_sv_bedpe -add_simple_gene_names", true );
		}
	}
	
	if ($use_dragen && ! $single_sample)
	{
		list($server) = exec2("hostname -f");
		//DRAGEN OUTFILES
		$dragen_output_vcf = "$dragen_output_folder/{$prefix}_dragen.vcf.gz";
		$dragen_output_msi = "$dragen_output_folder/{$prefix}_dragen_msi.json";
		$dragen_output_svs = "$dragen_output_folder/{$prefix}_dragen_svs.vcf.gz";
		$dragen_log_file = "$dragen_output_folder/{$prefix}_dragen.log";
		$sge_logfile = date("YmdHis")."_".implode("_", $server)."_".getmypid();
		$sge_update_interval = 300; //5min
		
		
		$t_bam_dragen = $t_bam;
		$n_bam_dragen = $n_bam;
		if (contains($t_bam, "/tmp/"))
		{
			#running based on a local bam file -> copy it to be reachable form dragen:
			$t_bam_dragen = $dragen_input_folder.basename($t_bam);
			$parser->copyFile($t_bam, $t_bam_dragen);
		}
		
		if (contains($n_bam, "/tmp/"))
		{
			#running based on a local bam file -> copy it to be reachable form dragen:
			$n_bam_dragen = $dragen_input_folder.basename($n_bam);
			$parser->copyFile($n_bam, $n_bam_dragen);
		}
		
		// create cmd for vc_dragen_somatic.php
		$args = array();
		$args[] = "-t_bam ".$t_bam_dragen;
		$args[] = "-out ".$dragen_output_vcf;
		// $args[] = "-out_sv ".$dragen_output_svs;
		$args[] = "-build ".$sys['build'];
		$args[] = "--log ".$dragen_log_file;
		
		if (! $single_sample)
		{
			$args[] = "-n_bam ".$n_bam_dragen;
		}
		
		if ($sys['type'] != "WGS")
		{
			$args[] = "-is_targeted";
		}
		
		$cmd = "php ".realpath(repository_basedir())."/src/NGS/vc_dragen_somatic.php ".implode(" ", $args);
		
		// submit GridEngine job to dragen queue
		$dragen_queues = explode(",", get_path("queues_dragen"));
		$sge_args = array();
		$sge_args[] = "-V";
		$sge_args[] = "-b y"; // treat as binary
		$sge_args[] = "-wd $dragen_output_folder";
		$sge_args[] = "-m n"; // switch off messages
		$sge_args[] = "-M ".get_path("queue_email");
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
			list($stdout) = exec2("qstat -u '*' | egrep '^\s+{$sge_id}\s+' 2>&1", false);
			$finished = trim(implode("", $stdout))=="";

			// log running state
			if (!$finished)
			{
				$state = explode(" ", preg_replace('/\s+/', ' ', $stdout[0]))[5];
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
		//copy vcf
		$parser->moveFile($dragen_output_vcf, $dragen_call_folder.basename($dragen_output_vcf));
		$parser->moveFile($dragen_output_vcf.".tbi", $dragen_call_folder.basename($dragen_output_vcf).".tbi");
		
		if (file_exists($dragen_output_msi))
		{
			$parser->moveFile($dragen_output_msi, $dragen_call_folder.basename($dragen_output_msi));
		}
		
		//filter dragen vcf
		$args = array();
		$args[] = "-in ".$dragen_call_folder.basename($dragen_output_vcf);
		$args[] = "-tumor_name {$t_id}";
		$args[] = "-normal_name {$n_id}";
		$args[] = "-out {$variants}";
		$args[] = "-build ".$sys['build'];
		$parser->execTool("NGS/an_filter_dragen_somatic.php", implode(" ", $args));
		
		
		if (is_file($dragen_output_svs))
		{
			//copy svs
			$parser->moveFile($dragen_output_svs, $dragen_call_folder.basename($dragen_output_svs));
			$parser->moveFile($dragen_output_svs.".tbi", $dragen_call_folder.basename($dragen_output_svs).".tbi");
				
			$parser->copyFile($dragen_call_folder.basename($dragen_output_svs), $manta_sv);
			$parser->copyFile($dragen_call_folder.basename($dragen_output_svs).".tbi", $manta_sv.".tbi");
				
			exec2(get_path("ngs-bits") . "VcfToBedpe -in $manta_sv -out $manta_sv_bedpe");
			if(!$single_sample)
			{
				$parser->execTool("Tools/bedpe2somatic.php", "-in $manta_sv_bedpe -out $manta_sv_bedpe -tid $t_id -nid $n_id");
			}
			
			if( db_is_enabled("NGSD") )
			{
				$parser->exec(get_path("ngs-bits") . "BedpeGeneAnnotation", "-in $manta_sv_bedpe -out $manta_sv_bedpe -add_simple_gene_names", true );
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
	else //NO dragen calling
	{
		if ($single_sample)
		{
			$parser->execTool("NGS/vc_varscan2.php", "-bam $t_bam -out $variants -build " .$sys['build']. " -target ". $roi. " -name $t_id -min_af $min_af");
			$parser->exec("tabix", "-f -p vcf $variants", true);
		}
		else
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
			$parser->execTool("NGS/vc_strelka2.php", implode(" ", $args_strelka));
		}
	}

	//add somatic BAF file
	if ($single_sample)
	{
		//create b-allele frequency file
		$params = array();
		$params[] = "-vcf {$variants}";
		$params[] = "-name {$t_id}";
		$params[] = "-out {$ballele}";
		
		if ($sys['type'] === "WGS")
		{
			$baf_args[] = "-downsample 100";
		}
		
		$parser->execTool("NGS/baf_germline.php", implode(" ", $params));
	}
	else
	{
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
			
			$parser->execTool("NGS/baf_somatic.php", implode(" ", $baf_args));
		}
		else
		{
			trigger_error("Cannot create BAF file as normal VCF missing: {$variants_germline_vcf}", E_USER_NOTICE);
		}
	}
	
	if (file_exists($variants) && !$skip_signatures && !$single_sample)
	{
		$snv_signatures_out = $out_folder."/snv_signatures/";
		$tmp_variants = $parser->tempFile(".vcf", "snv_signatures_");
		$parser->exec("bgzip","-c -d -@ {$threads} $variants > {$tmp_variants}", true);
		$parser->exec(get_path("python3")." ".repository_basedir()."/src/NGS/extract_signatures.py", "--in {$tmp_variants} --mode snv --outFolder {$snv_signatures_out} --reference GRCh38 --threads {$threads}", true);
	}
}

//Viral sequences alignment
$viral         = "{$t_basename}_viral.tsv";					// viral sequences results
$viral_bam     = "{$t_basename}_viral.bam";					// viral sequences alignment
$viral_bam_raw = "{$t_basename}_viral_before_dedup.bam";		// viral sequences alignment (no deduplication)
$viral_bed     = get_path("data_folder") . "/enrichment/somatic_viral.bed"; //viral enrichment
$viral_genome  = get_path("data_folder") . "/genomes/somatic_viral.fa"; //viral reference genome
if (in_array("vi", $steps))
{
	if(!file_exists($viral_genome) || !file_exists($viral_bed))
	{
		trigger_error("Could not find viral reference genome {$viral_genome} or target file {$viral_bed}. Skipping step \"vi\".", E_USER_WARNING);
	}
	else
	{
		//detection of viral sequences
		$t_bam_dedup = "{$t_basename}_before_dedup.bam";
		$t_bam_map_qc = "{$t_basename}_stats_map.qcML";
		$dedup_used = file_exists($t_bam_dedup);
		$vc_viral_args = [
			"-in ".($dedup_used ? $t_bam_dedup : $t_bam),
			"-viral_bam ".$viral_bam,
			"-viral_bam_raw ".$viral_bam_raw,
			"-viral_cov ".$viral,
			"-viral_chrs chrNC_007605",
			"-build_viral somatic_viral",
			"-in_qcml ".$t_bam_map_qc,
			"-threads ".$threads
		];
		if ($dedup_used)
		{
			$vc_viral_args[] = "-barcode_correction";
		}
		$parser->execTool("NGS/vc_viral_load.php", implode(" ", $vc_viral_args));
	}
}

//CNV calling
$som_clincnv = $full_prefix . "_clincnv.tsv"; //ClinCNV output file
if(in_array("cn",$steps))
{
	// copy number variant calling
	$tmp_folder = $parser->tempFolder();
	
	//Generate file with off-target region
	$off_target_bed = get_path("data_folder")."/coverage/off_target_beds/".$sys['name_short'].".bed";
	if(!file_exists($off_target_bed))	create_off_target_bed_file($off_target_bed,$sys['target_file'], $ref_genome);
	
	/***************************************************
	 * GENERATE AND COPY COVERAGE FILES TO DATA FOLDER *
	 ***************************************************/
	 
	//directory with reference coverage files
	$ref_folder_t = get_path("data_folder")."/coverage/".$sys['name_short']."-tumor";

	 
	// coverage for tumor sample
	$t_cov = "{$tmp_folder}/{$t_id}.cov";
	$ref_file_t = "{$ref_folder_t}/{$t_id}.cov";
	
	if(!file_exists($ref_file_t) )
	{
		$parser->exec(get_path("ngs-bits")."BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam $t_bam -in $roi -out $t_cov -threads {$threads} -ref {$ref_genome}",true);
		$parser->exec(get_path("ngs-bits")."BedSort", "-uniq -in $t_cov -out $t_cov",true);
	}
	else 
	{
		$t_cov = $ref_file_t;
	}
	
	//directory with reference off-target coverage files
	$ref_folder_t_off_target = $ref_folder_t . "_off_target";
	
	// coverage for tumor sample (off-target)
	$t_cov_off_target = "{$tmp_folder}/{$t_id}_off_target.cov";
	$ref_file_t_off_target = "{$ref_folder_t_off_target}/{$t_id}.cov";
	
	if( !file_exists($ref_file_t_off_target ) )
	{
		$parser->exec(get_path("ngs-bits")."BedCoverage", "-clear -min_mapq 10 -decimals 4 -in $off_target_bed -bam $t_bam -out $t_cov_off_target -threads {$threads} -ref {$ref_genome}",true);
		$parser->exec(get_path("ngs-bits")."BedSort", "-uniq -in $t_cov_off_target -out $t_cov_off_target",true);
	}
	else 
	{
		$t_cov_off_target = $ref_file_t_off_target;
	}
	
	//folders with tumor reference coverage files (of same processing system)


	//copy tumor sample coverage file to reference folder (has to be done before ClinCNV call to avoid analyzing the same sample twice)
	if (db_is_enabled("NGSD") && is_valid_ref_tumor_sample_for_cnv_analysis($t_id))
	{
		//create reference folder if it does not exist
		create_directory($ref_folder_t);
		
		//copy file
		if(!file_exists($ref_file_t)) //Do not overwrite existing reference files in cov folder
		{
			$parser->copyFile($t_cov, $ref_file_t); 
			$t_cov = $ref_file_t;
		}
		
		//create reference folder for tumor off target coverage files
		create_directory($ref_folder_t_off_target);
		
		if(!file_exists($ref_file_t_off_target))
		{
			$parser->copyFile($t_cov_off_target,$ref_file_t_off_target);
			$t_cov_off_target = $ref_file_t_off_target;
		}
	}

	if($single_sample) //use clincnv_germline with
	{		

		//get directory for normal and Offtarget coverages
		$ref_folder_n = get_path("data_folder")."/coverage/".$sys['name_short'];
		create_directory($ref_folder_n);

		$ref_folder_n_off_target = $ref_folder_n . "_off_target";
		create_directory($ref_folder_n_off_target);
		
		//create BED file with GC and gene annotations - if missing
		$bed = $ref_folder_t."/roi_annotated.bed";
		if (!file_exists($bed))
		{
			$ngsbits = get_path("ngs-bits");
			$pipeline = [
					["{$ngsbits}BedAnnotateGC", "-in ".$sys['target_file']." -clear -ref ".genome_fasta($sys['build'])],
					["{$ngsbits}BedAnnotateGenes", "-out {$bed}"],
				];
			$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
		}

		//create BAF folder
		$baf_folder = get_path("data_folder")."/coverage/". $sys['name_short']."_bafs";
		create_directory($baf_folder);
		$baf_file = "{$baf_folder}/{$t_id}.tsv";
		//create BAf file if not available
		$error = False;
		if(!file_exists($baf_file))
		{
			$t_gsvar = "{$t_basename}.GSvar";
			create_baf_file($t_gsvar, $t_bam, $baf_file, $ref_genome, $error);
		}

		//perform CNV analysis		
		$args = array(
			"-cov {$t_cov}",
			"-cov_folder {$ref_folder_n}",
			"-bed {$bed}",
			"-out {$som_clincnv}",
			"-tumor_only",
			"-max_cnvs 200",
			"-bed_off {$off_target_bed}",
			"-cov_off {$t_cov_off_target}",
			"-cov_folder_off {$ref_folder_n_off_target}"
		);
		if(!$error)
		{
			$args[] = "-baf_folder {$baf_folder}";
		}
		$args[] = "--log ".$parser->getLogFile();
		$parser->execTool("NGS/vc_clincnv_germline.php", implode(" ", $args), true);

		// annotate CNV file
		$repository_basedir = repository_basedir();
		$data_folder = get_path("data_folder");
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$som_clincnv}", true);
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed -no_duplicates -url_decode -out {$som_clincnv}", true);
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2024-02.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$som_clincnv}", true);

		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2023_3.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$som_clincnv}", true);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$omim_file} -no_duplicates -url_decode -out {$som_clincnv}", true);
		}
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$repository_basedir}/data/gene_lists/genes.bed -no_duplicates -url_decode -out {$som_clincnv}", true);

		//annotate additional gene info
		$parser->exec(get_path("ngs-bits")."CnvGeneAnnotation", "-in {$som_clincnv}  -add_simple_gene_names -out {$som_clincnv}", true);
		// skip annotation if no connection to the NGSD is possible
		if (db_is_enabled("NGSD"))
		{
			//annotate overlap with pathogenic CNVs
			$parser->exec(get_path("ngs-bits")."NGSDAnnotateCNV", "-in {$som_clincnv} -out {$som_clincnv}", true);
		}
	}
	else //ClinCNV for differential sample
	{
		//reference directory for normal coverage files
		$ref_folder_n = get_path("data_folder")."/coverage/".$n_sys['name_short'];
		
		// coverage for normal sample
		$ref_file_n = $ref_folder_n."/".$n_id.".cov";
		$n_cov = "{$tmp_folder}/{$n_id}.cov";
		
		if(!file_exists($ref_file_n)) 
		{
			$parser->exec(get_path("ngs-bits")."BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam $n_bam -in ".$n_sys['target_file']." -out $n_cov -threads {$threads} -ref {$ref_genome}", true);
			$parser->exec(get_path("ngs-bits")."BedSort", "-uniq -in $n_cov -out $n_cov",true);
		}
		else 
		{
			$n_cov = $ref_file_n;
		}
		
		//reference directory for off_target normal coverage files
		$ref_folder_n_off_target = $ref_folder_n . "_off_target";
		
		// coverage for normal sample (off-target)
		$ref_file_n_off_target = "{$ref_folder_n_off_target}/{$n_id}.cov";
		$n_cov_off_target = "{$tmp_folder}/{$n_id}_off_target.cov";
		
		if(!file_exists($ref_file_n_off_target) )
		{
			$parser->exec(get_path("ngs-bits")."BedCoverage", "-clear -min_mapq 10 -decimals 4 -in $off_target_bed -bam $n_bam -out $n_cov_off_target -threads {$threads} -ref {$ref_genome}",true);
			$parser->exec(get_path("ngs-bits")."BedSort", "-uniq -in $n_cov_off_target -out $n_cov_off_target",true);
		}
		else 
		{
			$n_cov_off_target = $ref_file_n_off_target;
		}
		
		// copy normal sample coverage file to reference folder (only if valid and not yet there).
		if (db_is_enabled("NGSD") && is_valid_ref_sample_for_cnv_analysis($n_id))
		{
			//create reference folder if it does not exist
			
			create_directory($ref_folder_n);

			//copy file
			if (!file_exists($ref_file_n)) //do not overwrite existing coverage files in ref folder
			{
				$parser->copyFile($n_cov, $ref_file_n);
				$n_cov = $ref_file_n;
			}
			
			
			create_directory($ref_folder_n_off_target);
			
			if(!file_exists($ref_file_n_off_target))
			{
				$parser->copyFile($n_cov_off_target,$ref_file_n_off_target);
				$n_cov_off_target = $ref_file_n_off_target;
			}
		}
		
		//append tumor-normal IDs to list with tumor normal IDs (stored in same folder as tumor coverage files)
		$t_n_list_file = $ref_folder_t . "/" . "list_tid-nid.csv";

		// create folder
		if (!file_exists($ref_folder_t))
		{
			mkdir($ref_folder_t);
			// check if successfull
			if (!file_exists($ref_folder_t)) trigger_error("Couldn't create folder '$ref_folder_t'!", E_USER_ERROR);
		}
		
		if (!file_exists($t_n_list_file))
		{
			$header = "##THIS FILE CONTAINS TUMOR AND NORMAL IDS OF PROCESSING SYSTEM ".$n_sys['name_short']."\n";
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
		
		//Append tumor-normal pair to csv list in tumor coverage folder if not included yet
		$already_in_list = false;
		foreach(file($t_n_list_file) as $line)
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
		if(!$already_in_list)
		{
			file_put_contents($t_n_list_file,"{$t_id},{$n_id}\n", FILE_APPEND | LOCK_EX);
		}

		$baf_folder = get_path("data_folder")."/coverage/". $sys['name_short']."_bafs";
		create_directory($baf_folder);
		if(db_is_enabled("NGSD"))
		{
			$db_conn = DB::getInstance("NGSD");
			//Create BAF file for each sample with the same processing system if not existing
			complement_baf_folder($t_n_list_file,$baf_folder,$db_conn, $ref_genome);
		}
		else
		{
			$normal_gsvar = "{$n_basename}.GSvar";
			create_baf_file($normal_gsvar,$n_bam,"{$baf_folder}/{$n_id}.tsv", $ref_genome);
			create_baf_file($normal_gsvar,$t_bam,"{$baf_folder}/{$t_id}.tsv", $ref_genome);
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
			"-bed", $sys['target_file'],
			"-t_cov_off", $t_cov_off_target,
			"-n_cov_off", $n_cov_off_target,
			"-bed_off", $off_target_bed,
			"-baf_folder", $baf_folder,
			"-cohort_folder", $cohort_folder,
			"-threads {$threads}"
			];
			
			if($sys['type'] == "WES")
			{
				$args_clincnv[] = "-lengthS 9";
				$args_clincnv[] = "-scoreS 200";
				$args_clincnv[] = "-filterStep 2";
			}
			
			
			if(isset($cnv_baseline_pos))
			{
				$args_clincnv[] = "-guide_baseline $cnv_baseline_pos";
			}
			
			$parser->execTool("NGS/vc_clincnv_somatic.php",implode(" ",$args_clincnv));
			
			//Annotate cytoband and data from network of cancer genes
			$parser->execTool("NGS/an_somatic_cnvs.php","-cnv_in $som_clincnv -out $som_clincnv -include_ncg -include_cytoband");
		}
		else
		{
			trigger_error("Not enough reference tumor-normal coverage files for processing system {$system} found. Skipping CNV calling.\n", E_USER_NOTICE);
		}
	}
	
	//calculate HRD based on clincnv.tsv file
	if (file_exists($som_clincnv))
	{
		if (!$single_sample && !$skip_HRD)
		{
			$parser->execTool("NGS/an_scarHRD.php" , "-cnvs {$som_clincnv} -tumor {$t_id} -normal {$n_id} -out_folder {$out_folder}");
		}
		
		if(!$single_sample && !$skip_signatures)
		{
			$cnv_signatures_out = $out_folder."/cnv_signatures/";
			$parser->exec(get_path("python3")." ".repository_basedir()."/src/NGS/extract_signatures.py", "--in {$som_clincnv} --mode cnv --outFolder {$cnv_signatures_out} --reference GRCh38 --threads {$threads}", true);
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
	$parser->execTool("Pipelines/annotate.php", "-out_name $prefix -out_folder $tmp_folder1 -system $system -vcf $variants -somatic -threads $threads");

	// run somatic QC
	if (!$single_sample)
	{
		$links = array_filter([
			"{$t_basename}_stats_fastq.qcML",
			"{$t_basename}_stats_map.qcML",
			"{$n_basename}_stats_fastq.qcML",
			"{$n_basename}_stats_map.qcML"
		], "file_exists");

		$args_somaticqc = [
			"-tumor_bam", $t_bam,
			"-normal_bam", $n_bam,
			"-somatic_vcf", $tmp_vcf,
			"-target_bed", $roi,
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
		}

		$parser->exec(get_path("ngs-bits")."SomaticQC", implode(" ", $args_somaticqc), true);
	}

	//add sample info to VCF header
	$s = Matrix::fromTSV($tmp_vcf);
	$comments = $s->getComments();
	$comments[] = gsvar_sample_header($t_id, array("IsTumor" => "yes"), "#", "");
	if (!$single_sample)
	{
		$comments[] = gsvar_sample_header($n_id, array("IsTumor" => "no"), "#", "");
	}
	$s->setComments(sort_vcf_comments($comments));
	$s->toTSV($tmp_vcf);

	// zip and index vcf file
	$parser->exec("bgzip", "-c $tmp_vcf > $variants_annotated", true);
	$parser->exec("tabix", "-f -p vcf $variants_annotated", true);

	// convert vcf to GSvar
	$args = array("-in $tmp_vcf", "-out $variants_gsvar", "-t_col $t_id");
	if (!$single_sample) $args[] = "-n_col $n_id";
	$parser->execTool("NGS/vcf2gsvar_somatic.php", implode(" ", $args));
	
	//Annotate data from network of cancer genes
	$parser->execTool("NGS/an_somatic_gsvar.php" , "-gsvar_in $variants_gsvar -out $variants_gsvar -include_ncg");

	//Determine cfDNA monitoring candidates (only tumor-normal samples)
	$umiVar2_path = get_path("umiVar2");
	if (!$single_sample)
	{
		if (file_exists($umiVar2_path."/select_monitoring_variants.py"))
		{
			//remove previous calls
			if (file_exists($cfdna_folder)) exec2("rm -r $cfdna_folder"); 
			//set parameters
			$params = array();
			$params[] = "-v ${tmp_vcf}";
			$params[] = "-g ${variants_gsvar}";
			$params[] = "-o ${cfdna_folder}";
			$params[] = "-r ${ref_genome}";
			$params[] = "-n 65"; // select 50 candidate variants
			$params[] = "-i"; // ignore INDELS
			// call variant selection in virtual environment
			//set environment variables
			putenv("umiVar_python_binary=\"".get_path("python3")."\"");
			putenv("umiVar_R_binary=\"".get_path("rscript")."\"");
			putenv("umiVar_samtools_binary=\"".get_path("samtools")."\"");
			$parser->exec(get_path("python3"), $umiVar2_path."/select_monitoring_variants.py ".implode(" ", $params));
		}
		else
		{
			trigger_error("UmiVar2 cannot be found! Cannot preselect variants for cfDNA analysis!", E_USER_WARNING);
		}
	}
	else
	{
		trigger_error("cfDNA candidate preselection only supported for tumor-normal samples.", E_USER_NOTICE);
	}
}

//MSI calling
$msi_o_file = $full_prefix . "_msi.tsv";						//MSI
if (in_array("msi", $steps) && !$single_sample)
{
	//file that contains MSI in target
	$msi_ref = get_path("data_folder") . "/dbs/msisensor-pro/msisensor_references_".$n_sys['build'].".list";
	
	if(!file_exists($msi_ref))
	{
		print("Could not find loci reference file $msi_ref. Trying to generate it.\n");
		$parser->exec(get_path("msisensor")," scan -d $ref_genome -o $msi_ref", false);
	}

	$parameters = "-n_bam $n_bam -t_bam $t_bam -msi_ref $msi_ref -threads $threads -out " .$msi_o_file. " -build ".$n_sys['build'];
	
	$parser->execTool("NGS/detect_msi.php",$parameters);
}
elseif ($single_sample && in_array("msi",$steps))
{
	trigger_error("Calling microsatellite instabilities is only possible for tumor normal pairs",E_USER_NOTICE);
}


if (in_array("an_rna", $steps))
{
	list($t_name) = explode("_", $t_id);

	//Determine tumor rna sample from NGSD via sample_relations
	$ps_rna_bams = array();

	if( db_is_enabled("NGSD") && !isset($t_rna_bam) )
	{
		$db = DB::getInstance("NGSD");
		$res = $db->executeQuery("SELECT id, name FROM sample WHERE name=:name", array("name" => $t_name));
		$sample_id = $res[0]["id"];
		
		//get related RNA samples from NGSD
		if (count($res) >= 1)
		{
			$relations = $db->executeQuery("SELECT * FROM sample_relations WHERE relation='same sample' AND (sample1_id=:sid OR sample2_id=:sid)", array("sid" => $sample_id));
			foreach($relations as $row)
			{
				$sample_id_annotation = $row['sample1_id'] != $sample_id ? $row['sample1_id'] : $row['sample2_id'];
				
				$res = $db->executeQuery("SELECT ps.sample_id, ps.process_id, ps.processing_system_id, ps.quality, sys.id, sys.type, CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as psample FROM processed_sample as ps, processing_system as sys, sample as s WHERE ps.sample_id=:sid AND sys.type='RNA' AND ps.processing_system_id=sys.id AND ps.sample_id=s.id", array("sid" => $sample_id_annotation));
				$rna_ids = array_column($res, 'psample');
				foreach($rna_ids as $rna_id)
				{
					$ps_rna_bams[$rna_id] = get_processed_sample_info($db, $rna_id)["ps_bam"];
				}
			}
		}
	}
	elseif(isset($t_rna_bam))
	{
		$ps_rna_bams[basename($t_rna_bam)] = $t_rna_bam; 
	}


	if(count($ps_rna_bams) < 1)
	{
		if (count($steps) > 1)
		{
			trigger_error("Skipping step an_rna!\nCouldn't find tumor RNA bam file. For annotation step \"an_rna\" tumor RNA bam file must be specified (via paramter -t_rna_bam or determined via sample_relations).", E_USER_WARNING);
			
			// remove an_rna step
			$key = array_search("an_rna", $steps);
			if ($key !== false) {
				unset($steps[$key]);
			}
		}
		else
		{
			trigger_error("Couldn't find tumor RNA bam file. For annotation step \"an_rna\" tumor RNA bam file must be specified (via paramter -t_rna_bam or determined via sample_relations).", E_USER_ERROR);
		}
	}
}

//RNA annotation
if (in_array("an_rna", $steps))
{	
	//Determine reference tissue type 1.) from parameter -rna_ref_tissue or 2.) from NGSD (if available)
	if(!isset($rna_ref_tissue) && !db_is_enabled("NGSD"))
	{
		trigger_error("For annotation step \"an_rna\" a tissue type for RNA reference data has to be specified.", E_USER_ERROR);
	}
	if(!isset($rna_ref_tissue) && db_is_enabled("NGSD"))
	{
		$db = DB::getInstance("NGSD");
		$tumor_db_id = get_processed_sample_id($db, $t_id);
		$s_id = $db->getValue("SELECT sample_id FROM processed_sample where id=$tumor_db_id");
		$res = $db->getValues("SELECT DISTINCT sdi.disease_info FROM sample s LEFT JOIN sample_relations sr ON s.id=sr.sample1_id OR s.id=sr.sample2_id LEFT JOIN sample_disease_info sdi ON sdi.sample_id=sr.sample1_id OR sdi.sample_id=sr.sample2_id WHERE s.id=$s_id AND sdi.type='RNA reference tissue' AND (sr.relation='same sample' OR sr.relation IS NULL)");
		if(count($res) == 1)
		{
			list($rna_ref_tissue) = $res;
		}
		else
		{
			trigger_error("Found multiple or no RNA reference tissue in NGSD. Aborting...", E_USER_ERROR);
		}
	}
	
	$rna_counts = array(); //file contains transcript counts
	foreach($ps_rna_bams as $rna_id => $rna_bam)
	{
		$rna_counts_tmp = glob(dirname($rna_bam)."/*_counts.tsv");
		if(count($rna_counts_tmp) != 1)
		{
			trigger_error("Could not find or found multiple RNA count files in sample folder " . dirname($rna_bam), E_USER_ERROR);
		}
		
		$rna_counts[$rna_id] = $rna_counts_tmp[0];
	}

	//Annotate data from all detected RNA files
	foreach($ps_rna_bams as $rna_id => $rna_bam)
	{
		check_genome_build($rna_bam, $sys['build']);
		$rna_count = $rna_counts[$rna_id];
		
		//Calculate sample similarity between tumor and RNA
		if ($skip_correlation)
		{
			trigger_error("The correlation calculation between DNA and RNA ({$rna_id}) is skipped!");
			continue;
		}
		$min_corr = 0.85;  //TODO: evaluate value
		
		if (file_exists($t_bam))
		{
			$output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", "-in {$rna_bam} {$t_bam} -mode bam -build ".ngsbits_build($sys['build']), true);
			$correlation = explode("\t", $output[0][1])[3];
			if ($correlation < $min_corr && ! $skip_correlation)
			{
				trigger_error("The genotype correlation of DNA and RNA ({$rna_id}) is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
			}
			else
			{
				trigger_error("The genotype correlation of DNA and RNA ({$rna_id}) is {$correlation}.", E_USER_NOTICE);
			}
		}
		else
		{
			trigger_error("BAM file does not exist for tumor DNA sample '{$t_bam}'!", E_USER_ERROR);
		}
		
		//SNVs
		$args = [
		"-gsvar_in $variants_gsvar",
		"-out $variants_gsvar",
		"-rna_id $rna_id",
		"-rna_counts $rna_count",
		"-rna_bam $rna_bam"];
		$parser->execTool("NGS/an_somatic_gsvar.php", implode(" ", $args));
		
		//CNVs
		if(file_exists($som_clincnv))
		{
			$parser->execTool("NGS/an_somatic_cnvs.php", " -cnv_in $som_clincnv -out $som_clincnv -rna_counts $rna_count -rna_id $rna_id -rna_ref_tissue " .str_replace(" ", 0, $rna_ref_tissue));
		}
	}
	
	//Reference tissue SNVs
	$args = [
		"-gsvar_in $variants_gsvar",
		"-out $variants_gsvar",
		"-rna_ref_tissue " .str_replace(" ", 0, $rna_ref_tissue)//Replace spaces by 0 because it is diffcult to pass spaces via command line.
	];
	$parser->execTool("NGS/an_somatic_gsvar.php", implode(" ", $args));

}

//Collect QC terms if necessary
$qc_other = $full_prefix."_stats_other.qcML";
if (in_array("vc", $steps) || in_array("vi", $steps) || in_array("msi", $steps) || in_array("cn", $steps) || in_array("db", $steps))
{
	$terms = array();
	$sources = array();
	
	// HRD score:
	$hrd_file = $full_prefix."_HRDresults.txt";
	if (file_exists($hrd_file))
	{
		foreach (file($hrd_file) as $line)
		{
			if (starts_with($line , '""')) continue;
			list($sample, $loh, $tai, $lst, $hrd) = explode("\t", trim($line));
			$terms[] = "QC:2000062\t{$loh}";
			$terms[] = "QC:2000063\t{$tai}";
			$terms[] = "QC:2000064\t{$lst}";
			$terms[] = "QC:2000126\t{$hrd}";
		}
		$sources[] = $hrd_file;
	}
	
	//virus:
	if (file_exists($viral))
	{
		$detected_viruses = [];
		foreach(file($viral) as $line)
		{
			if (starts_with($line, "#")) continue;
			
			list($chr, $start, $end, $v_name, $reads, $coverage, $coverage_rel, $mismatches, $ident) = explode("\t", $line);
			$v_base_name = explode("_", $v_name)[0];
			if (intval($reads) > 100 && ! in_array($v_base_name, $detected_viruses))
			{
				$detected_viruses[] = $v_base_name;
			}
		}
		
		$value = "None";
		if (count($detected_viruses) > 0)
		{
			$value = implode(", ", $detected_viruses);
		}
		
		$terms[] = "QC:2000130\t{$value}";
		$sources[] = $viral;
	}
	
	//MSI-status
	if (file_exists($msi_o_file))
	{
		foreach(file($msi_o_file) as $line)
		{
			if (trim($line) == "" || $line[0] == "#") continue;
			
			list($total, $somatic, $percent) = explode("\t", $line);
		}
		
		$terms[] = "QC:2000141\t{$percent}";
		$sources[] = $msi_o_file;
		
	}
	
	//HLA: TODO
	// if (file_exists($hla_file_tumor))
	// {
		
		// $sources[] = $hla_file_tumor;
	// }
	
	// if (!$single_sample && file_exists($hla_file_normal))
	// {
		
		// $sources[] = $hla_file_normal;
	// }
	
	//TODO mutational signatures
	
	//create qcML file
	if (count($sources) > 0) 
	{
		$tmp = $parser->tempFile("qc.tsv");
		file_put_contents($tmp, implode("\n", $terms));
		$parser->exec(get_path("ngs-bits")."TsvToQC", "-in $tmp -out $qc_other -sources ".implode(" ", $sources));
	}
}


//DB import
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	
	$t_info = get_processed_sample_info($db_conn, $t_id, false);
	$n_info = get_processed_sample_info($db_conn, $n_id, false);
	
	$ngsbits = get_path("ngs-bits");
	
	if (is_null($t_info) || (!$single_sample && is_null($n_info)))
	{
		trigger_error("No database import since no valid processing ID (T: {$t_id}".($single_sample ? "" : " /N: {$n_id}").")", E_USER_WARNING);
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
		$parser->exec("{$ngsbits}/NGSDImportSampleQC", "-ps $t_id -files $qcmls -force");

		// check tumor/normal flag
		if (!$t_info['is_tumor'])
		{
			trigger_error("Tumor sample $t_id is not flagged as tumor in NGSD!", E_USER_WARNING);
		}

		// analogous steps for normal sample, plus additional
		if (!$single_sample)
		{
			// check sex using control sample
			$parser->execTool("NGS/db_check_gender.php", "-in $n_bam -pid $n_id -sry_cov 30");

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

			
			// import SNVS (not for WGS)
			if (file_exists($variants_gsvar) && $sys['type'] !== "WGS")
			{
				check_genome_build($variants_gsvar, $sys['build']);
				$parser->exec(get_path("ngs-bits") . "/NGSDAddVariantsSomatic", " -t_ps $t_id -n_ps $n_id -var $variants_gsvar -var_force");
			}
			
			if(file_exists($som_clincnv) && $sys['type'] !== "WGS")
			{
				$parser->exec(get_path("ngs-bits") . "/NGSDAddVariantsSomatic", " -t_ps $t_id -n_ps $n_id -cnv $som_clincnv -cnv_force");
			}
			
			if(file_exists($manta_sv_bedpe) && $sys['type'] !== "WGS")
			{
				$parser->exec(get_path("ngs-bits") . "/NGSDAddVariantsSomatic", " -t_ps $t_id -n_ps $n_id -sv $manta_sv_bedpe -sv_force");
			}
		}
		
		//add secondary analysis (if missing)
		$parser->execTool("NGS/db_import_secondary_analysis.php", "-type 'somatic' -gsvar {$variants_gsvar}");
	}
}
?>
