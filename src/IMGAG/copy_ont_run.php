<?php 
/** 
	@page copy_ont_run
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$repo_folder = "/mnt/storage2/megSAP/pipeline"; //fixed absolute path to make the tests work for all users

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("copy_ont_run", "Copy ONT flow cell data (FASTQ, BAM).");
$parser->addString("run_dir",  "Input data folder.", false);

//optional
$parser->addString("run_name",  "Run name (will be derived from folder suffix if not provided.)", true, "");
$parser->addFlag("queue_sample", "Queue analysis of the sample.");
$parser->addFlag("queue_basecalling", "Queue basecalling of sample on GPU queue.");
$parser->addEnum("basecall_model",  "Model used for base calling: hac=high accuracy, sup=super accuracy", true, array("hac", "sup"), "hac");
$parser->addInt("min_qscore", "Minimal read qScore for basecalling", true, 9);
$parser->addFlag("force_basecalling", "Enforces complete (re-)basecalling of all POD5s in run folder.(requires '-queue_basecalling')");
$parser->addFlag("ignore_flowcell_id_check", "Ignores Errors due to mismatching flowcell IDs between NGSD and folder name (e.g. if flowcell was replaced during sequencing)");
$parser->addString("email",  "Email address which will be notified about the status of basecalling jobs", true, "");

//TODO: re-implement
// $parser->addFlag("fastq", "Copy FASTQ files.");
// $parser->addString("single_fastq",  "Create single FASTQ file, without sample lookup and database checks.", true, "");

$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addInt("threads", "Number of threads to use for file merging and compression.", true, 8);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));


/******** functions ********/

//check available data
function check_data_available($path, $subdir_name, $suffix)
{
	return file_exists($path) && is_dir($path) && (count(glob("{$path}/*/{$subdir_name}/*{$suffix}")) > 0);
}

//find subdirectory in directory
function get_sub_dirs($dir)
{
	$subdirs = array_values(array_diff(scandir($dir), array("..", ".")));
	if (count($subdirs) === 0)
	{
		trigger_error("No subdirectories found in '{$dir}'.", E_USER_ERROR);
	}
	foreach ($subdirs as &$subdir)
	{
		$subdir = "{$dir}/{$subdir}";
	}
	return $subdirs;
}

//check for skipped bases
function contains_skipped_bases($dirs, $threshold=5)
{
	$skipped_bases_found = false;
	foreach ($dirs as $dir)
	{
		if (is_dir("{$dir}/pod5_skip"))
		{
			list($stdout, $stderr, $exit_code) = exec2("du -b --summarize {$dir}/pod5_skip/");
			$folder_size = intval(explode("\t", $stdout[0])[0]);
			if ($folder_size > ($threshold*1024*1024)) // > $X MB
			{
				$skipped_bases_found = true;
			}
		}
	}

	return $skipped_bases_found;
}

//get run name
function get_run_name($dir)
{
	$tmp = explode("_", basename($dir));
	$suffix = end($tmp);
	if ((strlen($suffix) == 5) && preg_match('/^[0-9]+$/', $suffix))
	{
		return "#{$suffix}";
	}
	else
	{
		trigger_error("Invalid folder suffix '{$suffix}'! Cannot determine run name.", E_USER_ERROR);
	}
}

//get sample info from NGSD
function get_sample_info(&$db, $run_name)
{
	//find sample entered in database
	$result = $db->executeQuery(<<<SQL
	SELECT
		processed_sample.id,
		processed_sample.sequencing_run_id,
		(SELECT mid.name FROM mid WHERE mid.id = processed_sample.mid1_i7) as mid1_i7_name
	FROM
		processed_sample,
		sequencing_run
	WHERE
		sequencing_run.name = '{$run_name}' AND
		processed_sample.sequencing_run_id=sequencing_run.id
	SQL
	);

	if (count($result) < 1) trigger_error("No sample info found for run '{$run_name}'!", E_USER_ERROR);
 
	$sample_info_list = array();

	foreach ($result as $record)
	{
		$sample = processed_sample_name($db, $record["id"]);
		$sample_info = get_processed_sample_info($db, $sample);

		if (count($result) > 1)
		{
			preg_match('/_BP([0-9][0-9])Fw/', $record["mid1_i7_name"], $matches);
			if (count($matches) == 0)
			{
				trigger_error("Could not extract 2-digit barcode number from '{$record['mid1_i7_name']}'!", E_USER_ERROR);
			}
			$sample_info["barcode"] = "barcode{$matches[1]}";
		}

		$sample_info_list[$sample] = $sample_info;
	}
	return $sample_info_list;
}

//get BAM info
function get_bam_info($bam_file)
{
	$on_device_basecall_model = get_basecall_model($bam_file);

	$ret = execApptainer("samtools", "samtools view", "--count --exclude-flags 0x900 {$bam_file}", [$bam_file]);
	$num_records = intval($ret[0][0]);
	$ret = execApptainer("samtools", "samtools view", "--count --tag ML {$bam_file}", [$bam_file]);
	$num_base_mods = intval($ret[0][0]);
	$ret = execApptainer("samtools", "samtools view", "--count --require-flags 0x4 {$bam_file}", [$bam_file]);
	$num_unaligned = intval($ret[0][0]);

	$modified_bases = $num_records == $num_base_mods;
	$aligned = $num_unaligned < $num_records;

	return array($on_device_basecall_model, $modified_bases, $aligned);
}

//create basecall cmd
function create_basecall_command($run_dir, $out_bam, $basecall_model, $min_qscore, $secondary_output, $queue_sample, $skipped_only)
{
	$args = array();
	$args[] = "-run_dir {$run_dir}";
	$args[] = "-out_bam {$out_bam}";
	$args[] = "-basecall_model {$basecall_model}";
	$args[] = "-min_qscore {$min_qscore}";
	$args[] = "-secondary_output {$secondary_output}";

	if($queue_sample) $args[] = "-queue_sample";
	if($skipped_only) $args[] = "-skipped_only";
	$args[] = "-rename_skip_folder";
	$args[] = "--log ".dirname($run_dir)."/basecall_ont_run_".basename($run_dir)."_".date("YmdHis").".log";
	
	return "php ".realpath(repository_basedir())."/src/Tools/basecall_ont_run.php ".implode(" ", $args);
}

//create qsub cmd
function create_gpu_qsub_cmd($cmd_to_queue, $working_directory, $job_name, $email)
{
	list($server) = exec2("hostname -f");
	$sge_logfile = date("YmdHis")."_".implode("_", $server)."_".getmypid();

	$sge_args = array();
	$sge_args[] = "-V";
	$sge_args[] = "-b y"; // treat as binary
	$sge_args[] = "-wd {$working_directory}";
	if ($email == "") $sge_args[] = "-m n"; // switch off messages
	else 
	{
		$sge_args[] = "-m bea";
		$sge_args[] = "-M {$email}";
	}
	$sge_args[] = "-e ".get_path("gpu_log")."/$sge_logfile.err"; // stderr
	$sge_args[] = "-o ".get_path("gpu_log")."/$sge_logfile.out"; // stdout
	$sge_args[] = "-q ".get_path("queues_gpu"); // define queue
	$sge_args[] = "-N {$job_name}"; // set name
	$sge_args[] = "-pe smp 2"; //set slots
	$qsub_command_args = implode(" ", $sge_args)." ".$cmd_to_queue;

	return array("qsub", $qsub_command_args);

}

/******** init ********/

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile(dirname($run_dir)."/copy_ont_run_".basename($run_dir)."_".date("YmdHis").".log");


//absolute path
$run_dir = realpath($run_dir);
if (!file_exists($run_dir))
{
	trigger_error("Run directory '{$run_dir}' does not exist!", E_USER_ERROR);
}

//check parameter:
if ($force_basecalling && !$queue_basecalling) trigger_error("'-force_basecalling' requires '-queue_basecalling' to work!", E_USER_ERROR);

//set ulimit
$parser->exec("ulimit", "-n 10000");

//get sub directories
$subdirs = get_sub_dirs($run_dir);

// get run name
if ($run_name == "") $run_name = get_run_name($run_dir);

//database connection
$is_test_db = ($db == "NGSD_TEST");
$db = DB::getInstance($db);


/******** checks ********/

//check run name
$result = $db->executeQuery("SELECT id, fcid FROM sequencing_run WHERE name = '{$run_name}'");
if (count($result) !== 1)
{
	trigger_error("Sequencing run with the name {$run_name} not found in database.", E_USER_ERROR);
}
$run_id = $result[0]['id'];
$flowcell_id = $result[0]['fcid'];

//check flowcell ID
foreach ($subdirs as $subdir)
{
	if (!contains(basename($subdir), $flowcell_id))
	{
		trigger_error("Flowcell ID '{$flowcell_id}' not found in directory name '{$subdir}'.", ($ignore_flowcell_id_check)? E_USER_WARNING: E_USER_ERROR);
	}
}


/******** QC import ********/

//get sample info from NGSD
$sample_info = get_sample_info($db, $run_name);

//import run QC
$parser->execTool("Tools/runqc_parser_ont.php", "-name '{$run_name}' -run_dir {$run_dir} -force".(($is_test_db)?" -db NGSD_TEST":"").(($ignore_flowcell_id_check)?" -ignore_flowcell_id_check":""));

//GenLab import (skip on test run)
if (!$is_test_db)
{
	trigger_error("Importing information from GenLab...", E_USER_NOTICE);
	foreach ($sample_info as $sample => $info) 
	{
		$parser->execApptainer("ngs-bits", "NGSDImportGenlab", "-ps {$sample}");
	}
}




//TODO implement
// if ($fastq || $single_fastq)
if (false)
{
	//TODO: implement
	trigger_error("Not implemented yet!", E_USER_ERROR);

	

	// $fastq_available = check_data_available($run_dir, "fastq_pass", ".fastq.gz") || check_data_available($run_dir, "fastq_pass/barcode??", ".fastq.gz");
	// if (!$fastq_available && ($fastq || $single_fastq))
	// {
	// 	trigger_error("No FASTQ files available!", E_USER_ERROR);
	// }
	// $fastq_paths_glob = "{$run_dir}/*/fastq_pass";


	// if ($single_fastq !== "")
	// {
	// 	trigger_error("Copy and merge FASTQ files to given output file.", E_USER_NOTICE);

	// 	if (file_exists($single_fastq))
	// 	{
	// 		trigger_error("Output FASTQ file '{$single_fastq}' already exists, aborting.", E_USER_ERROR);
	// 	};

	// 	exec2("find {$fastq_paths_glob} -name '*.fastq.gz' -type f -exec zcat -f {} + | pigz -p {$threads} -c > {$single_fastq}");
	// 	trigger_error("FASTQ saved in {$single_fastq}.", E_USER_NOTICE);
	// 	exit();
	// }

	// if (!$bam_available && !$fastq_available && $prefer_bam)
	// {
	// 	trigger_error("No FASTQ/BAM files are available!", E_USER_ERROR);
	// }
	// $bam_paths_glob = "{$run_dir}/*/bam_pass";

	// if ($fastq || ($prefer_bam && !$bam_available))
	// {
	// 	trigger_error("Copy and merge FASTQ files.", E_USER_NOTICE);

	// 	$out_fastq = "{$out_dir}{$sample}.fastq.gz";
	// 	if (file_exists($out_fastq))
	// 	{
	// 		trigger_error("Output FASTQ file '{$out_fastq}' already exists, aborting.", E_USER_ERROR);
	// 	};

	// 	exec2("find {$fastq_paths_glob} -name '*.fastq.gz' -type f -exec cat  {} + > {$out_fastq}");
	// 	trigger_error("FASTQ saved in {$out_fastq}.", E_USER_NOTICE);	
	// }

	

}
else //bam output
{
	
	//check for pod5_skip subdirectory (partially basecalled flow cell)
	$contains_skipped_bases = contains_skipped_bases($subdirs);
	if ($contains_skipped_bases && !$queue_basecalling) trigger_error("'pod5_skip' directory present in '{$subdir}', some data has not been basecalled!", E_USER_ERROR);

	//check for barcodes
	$contains_barcodes = check_data_available($run_dir, "bam_pass/barcode??", ".bam");
	if ($contains_barcodes && $queue_basecalling) trigger_error("Automatic base calling with barcodes is currently not supported!", E_USER_ERROR);
	if ((count($sample_info) > 1) && $queue_basecalling) trigger_error("Automatic base calling with multiple samples per run is currently not supported!", E_USER_ERROR);


	$bam_available = check_data_available($run_dir, "bam_pass", ".bam") || check_data_available($run_dir, "bam_pass/barcode??", ".bam");
	if (!$bam_available && !$queue_basecalling) trigger_error("No BAM files available!", E_USER_ERROR);

	trigger_error("Run folder checks done", E_USER_NOTICE);

	foreach ($sample_info as $sample => $info) 
	{
		$bam_paths_glob = "{$run_dir}/*/bam_pass/*.bam";
		if (isset($info["barcode"]) && ($info["barcode"] != "")) $bam_paths_glob = "{$run_dir}/*/bam_pass/".$info["barcode"]."/*.bam";
		

		//get BAM info 
		

		$on_device_basecall_model = [];
		$modified_bases = [];
		$aligned = [];

		//check each subdir seperately
		foreach ($subdirs as $subdir) 
		{
			//get BAM files
			$bam_files = glob("{$subdir}/bam_pass/*.bam");
			if (count($bam_files) !== 0)
			{
				//check each subfolder
				list($on_device_basecall_model_subfolder, $modified_bases_subfolder, $aligned_subfolder) = get_bam_info($bam_files[array_rand($bam_files, 1)]);
				$on_device_basecall_model[$on_device_basecall_model_subfolder] = true;
				$modified_bases[$modified_bases_subfolder] = true;
				$aligned[$aligned_subfolder] = true;
			}
			else
			{
				// no bam files -> undefined basecall model
				$on_device_basecall_model["undefined"] = true;
				$modified_bases[false] = true;
				$aligned[false] = true;
			}
		}
		//special case: no basecalling
		if (count(array_keys($on_device_basecall_model)) == 0)
		{
			trigger_error("No BAMs found in subfolders!", E_USER_NOTICE);
			$on_device_basecall_model = "undefined";
			$modified_bases = false;
			$aligned = false;
		}
		else if (count(array_keys($on_device_basecall_model)) > 1)
		{
			trigger_error("Multiple basecall models found (".implode(", ", array_keys($on_device_basecall_model)).")! Cannot perform copy.", E_USER_ERROR);
		}
		else if (count(array_keys($modified_bases)) > 1)
		{
			trigger_error("Mixture of modified and unmodified bases found! Cannot perform copy.", E_USER_ERROR);
		}
		else
		{
			// default case: setting of all subfolders are identical
			$on_device_basecall_model = array_keys($on_device_basecall_model)[0];
			$modified_bases = array_keys($modified_bases)[0];
			$aligned = array_keys($aligned)[0];

			if ($aligned) trigger_error("Aligned BAM files available, but not supported yet. Will be treated as unaligned!", E_USER_WARNING);
			if ($modified_bases) trigger_error("Modified bases BAM files available.", E_USER_NOTICE);
		}

		//check on-device basecall model:
		if (str_contains($on_device_basecall_model, $basecall_model) && !$force_basecalling)
		{
			trigger_error("Basecall Models match: On-device basecall model: '{$on_device_basecall_model}', requested basecall model: '{$basecall_model}", E_USER_NOTICE);

			$out_dir = $info["ps_folder"];
				
			// copy to run folder during test:
			if ($is_test_db) $out_dir = $run_dir."/TEST_Sample_".$sample."/";
			
			//create sample folder 
			if (!file_exists($out_dir)) 
			{
				mkdir($out_dir, 0777, true);
				//apply file access permissions
				$parser->exec("chmod", "-R 775 {$out_dir}");
				$parser->exec("chgrp", "-R f_ad_bi_l_medgen_access_storages {$out_dir}", true, false);
			}
			$out_bam = "{$out_dir}{$sample}.unmapped.bam";
			if ($modified_bases) $out_bam = "{$out_dir}{$sample}.mod.unmapped.bam";
			if (file_exists($out_bam)) trigger_error("Output BAM file '{$out_bam}' already exists, aborting.", E_USER_ERROR);
			

			if ($contains_skipped_bases && $queue_basecalling)
			{

				//queue partial basecalling on GPU server
				trigger_error("Sample contains un-called bases. Queuing partial basecalling!", E_USER_NOTICE);

				$basecalling_cmd = create_basecall_command($run_dir, $out_bam, $basecall_model, $min_qscore, end($subdirs)."/bam_pass/", $queue_sample, true);
				$parser->log("Basecalling command:\t".$basecalling_cmd);

				$job_name = "megSAP_partial_basecalling_".basename($run_dir);
				$qsub_cmd = create_gpu_qsub_cmd($basecalling_cmd, $run_dir, $job_name, $email);

				// log sge command
				trigger_error("SGE command:\tqsub ".$qsub_cmd[1]);

				if ($is_test_db)
				{
					//only report qsub command
					trigger_error("Test mode, only reporting qsub command!", E_USER_WARNING);
				}
				else
				{
					// run qsub 
					list($stdout, $stderr) = $parser->exec($qsub_cmd[0], $qsub_cmd[1]);
					$sge_id = explode(" ", $stdout[0])[2];

					trigger_error("Basecalling queued with SGE id {$sge_id} and job name '{$job_name}'.", E_USER_NOTICE);
				}


			}
			else if (!$contains_skipped_bases)
			{
				//merge and copy to sample folder
				// concat bams
				$pipeline = [];
				$pipeline[] = ["find", "{$bam_paths_glob} -name '*.bam' -type f"];
				$pipeline[] = ["", $parser->execApptainer("samtools", "samtools cat", "-o {$out_bam} -b -", [$run_dir], [$out_dir], true)]; //no reference required
				$parser->execPipeline($pipeline, "merge unaligned BAM files");

				//apply file access permissions
				$parser->exec("chmod", "-R 775 {$out_dir}");
				$parser->exec("chgrp", "-R f_ad_bi_l_medgen_access_storages {$out_dir}", true, false);
				
				//queue sample
				if ($queue_sample) $parser->execTool("Tools/db_queue_analysis.php", "-samples {$sample} -type 'single sample'".(($is_test_db)?" -db NGSD_TEST -user unknown":""));
			}
			else trigger_error("Should not happen!", E_USER_ERROR);

		}
		else if ($queue_basecalling)
		{
			if ($force_basecalling) trigger_error("Full basecalling requested!", E_USER_NOTICE);
			else trigger_error("Basecall Models doesn't match: On-device basecall model: '{$on_device_basecall_model}', requested basecall model: '{$basecall_model}'. Full (re-)basecalling is needed!", E_USER_WARNING);

			//move basecalled files
			$basecall_suffix = "";
			if (str_contains($on_device_basecall_model, "hac")) $basecall_suffix = "_hac";
			else if (str_contains($on_device_basecall_model, "sup")) $basecall_suffix = "_sup";
			else if ($on_device_basecall_model == "undefined") $basecall_suffix = "_old";
			else trigger_error("Unknown basecall model type '{$on_device_basecall_model}'!", E_USER_ERROR);

			foreach ($subdirs as $subdir) 
			{
				if (is_dir($subdir."/bam_pass")) $parser->moveFile($subdir."/bam_pass", $subdir."/bam_pass".$basecall_suffix);
				if (is_dir($subdir."/bam_failed")) $parser->moveFile($subdir."/bam_failed", $subdir."/bam_failed".$basecall_suffix);
			}

			//create output folder (sample folder)
			$out_dir = $info["ps_folder"];
			// copy to run folder during test:
			if ($is_test_db) $out_dir = $run_dir."/TEST_Sample_".$sample."/";
			//create sample folder 
			if (!file_exists($out_dir)) 
			{
				mkdir($out_dir, 0777, true);
				//apply file access permissions
				$parser->exec("chmod", "-R 775 {$out_dir}");
				$parser->exec("chgrp", "-R f_ad_bi_l_medgen_access_storages {$out_dir}", true, false);
			}
			$out_bam = "{$out_dir}{$sample}.mod.unmapped.bam";
			// $out_bam = "{$out_dir}{$sample}.unmapped.bam";
			// if ($modified_bases) $out_bam = "{$out_dir}{$sample}.mod.unmapped.bam";
			if (file_exists($out_bam)) trigger_error("Output BAM file '{$out_bam}' already exists, aborting.", E_USER_ERROR);
			

			//queue full re-basecalling
			$basecalling_cmd = create_basecall_command($run_dir, $out_bam, $basecall_model, $min_qscore, end($subdirs)."/bam_pass/", $queue_sample, false);
			$parser->log("Basecalling command:\t".$basecalling_cmd);

			$job_name = "megSAP_basecalling_".basename($run_dir);
			$qsub_cmd = create_gpu_qsub_cmd($basecalling_cmd, $run_dir, $job_name, $email);

			if ($is_test_db)
			{
				//only report qsub command
				print "qsub command:\n".$qsub_cmd[0]." ".$qsub_cmd[1]."\n";
			}
			else
			{
				// log sge command
				$parser->log("SGE command:\tqsub ".$qsub_cmd[1]);

				// run qsub 
				list($stdout, $stderr) = $parser->exec($qsub_cmd[0], $qsub_cmd[1]);
				$sge_id = explode(" ", $stdout[0])[2];

				trigger_error("Basecalling queued with SGE id {$sge_id} and job name '{$job_name}'.", E_USER_NOTICE);
			}

		}
		else
		{
			trigger_error("Basecall Models doesn't match: On-device basecall model: '{$on_device_basecall_model}', requested basecall model: '{$basecall_model}'. Please check setting or re-do basecalling!", E_USER_ERROR);
		}

		
	}

}

if ($queue_sample)
{
	// update sequencing run analysis status
	$db->executeStmt("UPDATE sequencing_run SET sequencing_run.status='analysis_started' WHERE sequencing_run.name = '{$run_name}' AND sequencing_run.status IN ('run_started', 'run_finished')");
}

/*
	
	//check random BAM file to find out which information is present
	$bam_files = glob("{$bam_paths_glob}/*.bam");
	if (count($bam_files) !== 0)
	{
		$bam_random_file = $bam_files[array_rand($bam_files, 1)];
	
		$on_device_basecall_model = get_basecall_model($bam_random_file);
		trigger_error("Model used for on-device basecalling: {$on_device_basecall_model}", E_USER_NOTICE);

		$ret = $parser->execApptainer("samtools", "samtools view", "--count --exclude-flags 0x900 {$bam_random_file}", [$run_dir]);
		$num_records = intval($ret[0][0]);
		$ret = $parser->execApptainer("samtools", "samtools view", "--count --tag ML {$bam_random_file}", [$run_dir]);
		$num_base_mods = intval($ret[0][0]);
		$ret = $parser->execApptainer("samtools", "samtools view", "--count --require-flags 0x4 {$bam_random_file}", [$run_dir]);
		$num_unaligned = intval($ret[0][0]);

		$modified_bases = $num_records == $num_base_mods;
		$aligned = $num_unaligned < $num_records;

		if ($aligned) trigger_error("Aligned BAM files available.", E_USER_NOTICE);
		if ($modified_bases) trigger_error("Modified bases BAM files available.", E_USER_NOTICE);

		
	}
	else
	{
		$aligned = false;
		$modified_bases = false;
	}

}

foreach ($result as $record)
{
	$sample = processed_sample_name($db_con, $record["id"]);
	$sample_info = get_processed_sample_info($db_con, $sample);

	if (count($result) >= 2)
	{
		preg_match('/_BP([0-9][0-9])Fw/', $record["mid1_i7_name"], $matches);
		if (count($matches) == 0)
		{
			trigger_error("Could not extract 2-digit barcode number from '{$record['mid1_i7_name']}'!", E_USER_ERROR);
		}
		$barcode = "barcode{$matches[1]}";
		$bam_paths_glob = "{$run_dir}/* /bam_pass/{$barcode}";
		$fastq_paths_glob = "{$run_dir}/* /fastq_pass/{$barcode}";
	}

	$out_dir = $sample_info["ps_folder"];
	// copy to run folder during test:
	if ($db=="NGSD_TEST")
	{
		$out_dir = $run_dir."/TEST_Sample_".$sample."/";
	}

	if (!file_exists($out_dir))
	{
		mkdir($out_dir, 0777, true);
	}




	if (($bam_available && $prefer_bam) || $bam)
	{
		trigger_error("Copy and merge BAM files.", E_USER_NOTICE);
		$genome = genome_fasta($build);


		if ($aligned && !$ignore_aligned)
		{
			$out_bam = "{$out_dir}{$sample}.bam";
		}
		elseif ($modified_bases)
		{
			$out_bam = "{$out_dir}{$sample}.mod.unmapped.bam";
		}
		else
		{
			$out_bam = "{$out_dir}{$sample}.unmapped.bam";
		}

		if (file_exists($out_bam))
		{
			trigger_error("Output BAM file '{$out_bam}' already exists, aborting.", E_USER_ERROR);
		}

		$pipeline = [];
		$pipeline[] = ["find", "{$bam_paths_glob} -name '*.bam' -type f"];

		if ($aligned && !$ignore_aligned)
		{
			$tmp_for_sorting = $parser->tempFile();
			//merge presorted files
			$pipeline[] = ["", $parser->execApptainer("samtools", "samtools merge", "--reference {$genome} --threads {$threads} -b - -o {$out_bam}", [$genome, $run_dir], [$out_dir], true)];
			$parser->execPipeline($pipeline, "merge aligned BAM files");
			$parser->indexBam($out_bam, $threads);

		}
		else
		{
			$pipeline[] = ["", $parser->execApptainer("samtools", "samtools cat", "-o {$out_bam} -b -", [$run_dir], [$out_dir], true)]; //no reference required
			$parser->execPipeline($pipeline, "merge unaligned BAM files");
		}
	}

	
	//apply file access permissions
	exec2("chmod -R 775 {$out_dir}");
	exec2("chgrp -R f_ad_bi_l_medgen_access_storages {$out_dir}", false);

	if ($queue_sample)
	{
		$parser->execTool("Tools/db_queue_analysis.php", "-samples {$sample} -type 'single sample' -db {$db}".(($db=="NGSD_TEST")?" -user unknown":""));
	}
}

*/

?>