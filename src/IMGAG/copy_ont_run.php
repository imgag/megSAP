<?php
/**
 * @page copy_ont_run
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
//TODO: remove
//$repo_folder = "/mnt/storage2/megSAP/pipeline"; //fixed absolute path to make the tests work for all users

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("copy_ont_run", "Copy ONT flow cell data (BAM).");
$parser->addString("run_dir", "Input data folder.", false);

// optional
$parser->addString("run_name", "Run name (will be derived from folder suffix if not provided.)", true, "");
$parser->addFlag("queue_sample", "Queue analysis of the sample.");
$parser->addFlag("queue_basecalling", "Queue basecalling of sample on GPU queue.");
$parser->addEnum("basecall_model", "Model used for base calling: hac=high accuracy, sup=super accuracy", true, array("hac", "sup"), "hac");
$parser->addInt("min_qscore", "Minimal read qScore for basecalling", true, 9);
$parser->addFlag("force_basecalling", "Enforces complete (re-)basecalling of all POD5s in run folder.(requires '-queue_basecalling')");
$parser->addFlag("ignore_flowcell_id_check", "Ignores Errors due to mismatching flowcell IDs between NGSD and folder name (e.g. if flowcell was replaced during sequencing)");
$parser->addString("email", "Email address which will be notified about the status of basecalling jobs", true, "");
$parser->addFlag("basecalling_copy_pod5", "In case of basecalling: Copy pod5 files to /tmp to reduce IO/network load.");
$parser->addFlag("ignore_copy_complete_check", "Ignores Errors due to missing 'IMGAG_CopyDone.txt' file in run folder.");

$parser->addEnum("db", "Database to connect to.", true, db_names(), "NGSD");
$parser->addInt("threads", "Number of threads to use for file merging and compression.", true, 8);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));


/******** functions ********/

//check available data
function check_data_available($path, $subdir_name, $suffix)
{
	return is_dir($path) && (count(glob("{$path}/*/{$subdir_name}/*{$suffix}")) > 0);
}

//find subdirectory in directory
function get_sub_dirs($dir)
{
	$entries = array_values(array_diff(scandir($dir), array("..", ".", "IMGAG_CopyDone.txt")));
	$subdirs = array();

	foreach ($entries as $entry)
	{
		// ignore folders with leading '_':
		if (str_starts_with($entry, "_")) 
		{
			trigger_error("Skipping folder {$entry}.", E_USER_WARNING);
			continue;
		}

		$full = "{$dir}/{$entry}";
		if (is_dir($full))
		{
			$subdirs[] = $full;
		}
	}

	if (count($subdirs) === 0)
	{
		trigger_error("No subdirectories found in '{$dir}'.", E_USER_ERROR);
	}

	return $subdirs;
}

function contains_skipped_bases($dirs, $threshold = 5)
{
	foreach ($dirs as $dir)
	{
		$skip_dir = "{$dir}/pod5_skip";
		if (is_dir($skip_dir))
		{
			list($stdout, $stderr, $exit_code) = exec2("du -b --summarize {$skip_dir}");
			$folder_size = intval(explode("\t", $stdout[0])[0]);
			if ($folder_size > ($threshold * 1024 * 1024)) return true;
		}
	}

	return false;
}

//get run name
function get_run_name($dir)
{
	$tmp = explode("_", basename($dir));
	$suffix = end($tmp);

	if ((strlen($suffix) === 5) && preg_match('/^[0-9]+$/', $suffix))
	{
		return "#{$suffix}";
	}

	trigger_error("Invalid folder suffix '{$suffix}'! Cannot determine run name.", E_USER_ERROR);
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
			if (count($matches) === 0)
			{
				trigger_error("Could not extract 2-digit barcode number from '{$record['mid1_i7_name']}'!", E_USER_ERROR);
			}
			$sample_info["barcode"] = "barcode{$matches[1]}";
		}

		$sample_info_list[$sample] = $sample_info;
	}

	return $sample_info_list;
}

//create basecall cmd
function create_basecall_command($run_dir, $out_bam, $basecall_model, $min_qscore, $secondary_output, $queue_sample, $skipped_only, $basecalling_copy_pod5)
{
	$args = array();
	$args[] = "-run_dir {$run_dir}";
	$args[] = "-out_bam {$out_bam}";
	$args[] = "-basecall_model {$basecall_model}";
	$args[] = "-min_qscore {$min_qscore}";
	$args[] = "-secondary_output {$secondary_output}";

	if ($queue_sample) $args[] = "-queue_sample";
	if ($skipped_only) $args[] = "-skipped_only";
	if ($basecalling_copy_pod5) $args[] = "-copy_pod5";
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
	if ($email === "") $sge_args[] = "-m n"; // switch off messages
	else 
	{
		$sge_args[] = "-m bea";
		$sge_args[] = "-M {$email}";
	}

	$sge_args[] = "-e ".get_path("gpu_log")."/{$sge_logfile}.err";
	$sge_args[] = "-o ".get_path("gpu_log")."/{$sge_logfile}.out";
	$sge_args[] = "-q ".get_path("queues_gpu");
	$sge_args[] = "-N ".$job_name;
	$sge_args[] = "-pe smp 2";

	return array("qsub", implode(" ", $sge_args)." ".$cmd_to_queue);
}

function get_bam_files_for_sample($subdirs, $barcode = "")
{
	$bam_files = array();

	foreach ($subdirs as $subdir)
	{
		$pattern = "{$subdir}/bam_pass/*.bam";
		if ($barcode !== "")
		{
			$pattern = "{$subdir}/bam_pass/{$barcode}/*.bam";
		}

		foreach (glob($pattern) as $file)
		{
			if (is_file($file)) $bam_files[] = $file;
		}
	}

	sort($bam_files);
	return $bam_files;
}

/******** init ********/

//absolute path
$resolved_run_dir = realpath($run_dir);
if ($resolved_run_dir === false || !is_dir($resolved_run_dir))
{
	trigger_error("Run directory '{$run_dir}' does not exist!", E_USER_ERROR);
}
$run_dir = $resolved_run_dir;

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") 
{
	$parser->setLogFile(dirname($run_dir)."/copy_ont_run_".basename($run_dir)."_".date("YmdHis").".log");
}


//check parameter:
if ($force_basecalling && !$queue_basecalling) 
{
	trigger_error("'-force_basecalling' requires '-queue_basecalling' to work!", E_USER_ERROR);
}

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

//check if copy is complete
if (!file_exists($run_dir."/IMGAG_CopyDone.txt"))
{
	trigger_error("'IMGAG_CopyDone.txt' file missing! Copy from sequencer may not be completed.", ($ignore_copy_complete_check ? E_USER_WARNING : E_USER_ERROR));
}

$flowcell_id = $result[0]['fcid'];

//check flowcell ID
foreach ($subdirs as $subdir)
{
	if (!contains(basename($subdir), $flowcell_id))
	{
		trigger_error("Flowcell ID '{$flowcell_id}' not found in directory name '{$subdir}'.", ($ignore_flowcell_id_check ? E_USER_WARNING : E_USER_ERROR));
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


//check for pod5_skip subdirectory (partially basecalled flow cell)
$contains_skipped_bases = contains_skipped_bases($subdirs);
if ($contains_skipped_bases && !$queue_basecalling) trigger_error("'pod5_skip' directory present in '".implode(",", $subdirs)."', some data has not been basecalled!", E_USER_ERROR);

//check for barcodes
$contains_barcodes = check_data_available($run_dir, "bam_pass/barcode??", ".bam");
if ($contains_barcodes && $queue_basecalling) trigger_error("Automatic base calling with barcodes is currently not supported!", E_USER_ERROR);
if ((count($sample_info) > 1) && $queue_basecalling) trigger_error("Automatic base calling with multiple samples per run is currently not supported!", E_USER_ERROR);


$bam_available = check_data_available($run_dir, "bam_pass", ".bam") || check_data_available($run_dir, "bam_pass/barcode??", ".bam");
if (!$bam_available && !$queue_basecalling) trigger_error("No BAM files available!", E_USER_ERROR);

trigger_error("Run folder checks done", E_USER_NOTICE);

foreach ($sample_info as $sample => $info)
{
	$barcode = (isset($info["barcode"]) && $info["barcode"] != "") ? $info["barcode"] : "";
	$all_bam_files = get_bam_files_for_sample($subdirs, $barcode);

	//get BAM info 
	

	$on_device_basecall_model = [];
	$modified_bases = [];
	$aligned = [];

	//check each subdir seperately
	foreach ($subdirs as $subdir)
	{
		$bam_files = get_bam_files_for_sample(array($subdir), $barcode);

		if (count($bam_files) !== 0)
		{
			//check each subfolder
			$rng_bam_file = $bam_files[array_rand($bam_files, 1)];
			$model = get_basecall_model($rng_bam_file);
			$has_mods = contains_methylation($rng_bam_file, ($is_test_db ? 12 : 100), $build);
			$is_aligned = is_bam_aligned($rng_bam_file, $build);

			$on_device_basecall_model[$model] = true;
			$modified_bases[(int)$has_mods] = true;
			$aligned[(int)$is_aligned] = true;
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
	if ((count(array_keys($on_device_basecall_model)) === 1) && (array_keys($on_device_basecall_model)[0] === "undefined"))
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
		$on_device_basecall_model = array_keys($on_device_basecall_model)[0];
		$modified_bases = ((int)array_keys($modified_bases)[0] === 1);
		$aligned = ((int)array_keys($aligned)[0] === 1);

		if ($aligned) trigger_error("Aligned BAM files available, but not supported yet. Will be treated as unaligned!", E_USER_WARNING);
		if ($modified_bases) trigger_error("Modified bases BAM files available.", E_USER_NOTICE);
	}

	//check on-device basecall model:
	if (str_contains($on_device_basecall_model, $basecall_model) && !$force_basecalling)
	{
		trigger_error("Basecall models match: On-device basecall model: '{$on_device_basecall_model}', requested basecall model: '{$basecall_model}'", E_USER_NOTICE);

		$out_dir = $info["ps_folder"];
			
		// copy to run folder during test:
		if ($is_test_db) $out_dir = $run_dir."/TEST_Sample_".$sample."/";
		
		//create sample folder 
		if (!file_exists($out_dir)) 
		{
			mkdir($out_dir, 0775, true);
			//apply file access permissions
			$parser->exec("chgrp", "-R f_ad_bi_l_medgen_access_storages {$out_dir}", true, false);
		}
		if ($modified_bases) $out_bam = "{$out_dir}{$sample}.mod.unmapped.bam";
		else $out_bam = "{$out_dir}{$sample}.unmapped.bam";
		if (file_exists($out_bam)) trigger_error("Output BAM file '{$out_bam}' already exists, aborting.", E_USER_ERROR);
		

		if ($contains_skipped_bases && $queue_basecalling)
		{

			//queue partial basecalling on GPU server
			trigger_error("Sample contains un-called bases. Queuing partial basecalling!", E_USER_NOTICE);

			$basecalling_cmd = create_basecall_command($run_dir, $out_bam, $basecall_model, $min_qscore, end($subdirs)."/bam_pass/", $queue_sample, true, $basecalling_copy_pod5);
			$parser->log("Basecalling command:\t".$basecalling_cmd);

			$job_name = "megSAP_partial_basecalling_".basename($run_dir);
			$qsub_cmd = create_gpu_qsub_cmd($basecalling_cmd, $run_dir, $job_name, $email);

			// log sge command
			trigger_error("SGE command:\tqsub ".$qsub_cmd[1], E_USER_NOTICE);

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
			// merge and copy to sample folder
			if (count($all_bam_files) === 0)
			{
				trigger_error("No BAM files found for sample '{$sample}'.", E_USER_ERROR);
			}
			// create temp bam file list
			$tmp_bam_file_list = $parser->tempFile("bam_list.txt");
			file_put_contents($tmp_bam_file_list, implode("\n", $all_bam_files)."\n");

			// concat bams
			$parser->execApptainer("samtools", "samtools cat", "-o {$out_bam} -b {$tmp_bam_file_list}", [$run_dir], [$out_dir]); //no reference required

			//apply file access permissions
			$parser->exec("chmod", "-R 775 {$out_dir}");
			$parser->exec("chgrp", "-R f_ad_bi_l_medgen_access_storages {$out_dir}", true, false);
			
			//queue sample
			if ($queue_sample) 
			{
				$args = "-samples {$sample} -type 'single sample'";
				if ($is_test_db) $args .= " -db NGSD_TEST -user unknown";
				$parser->execTool("Tools/db_queue_analysis.php", $args);
				//check if trio needs to be queued
				$parents = get_trio_parents($db, $sample);
				if (!is_null($parents))
				{
					list($father, $mother) = $parents;
					$trio_parameters = "-type 'trio' -samples {$sample} {$father} {$mother} -info child father mother";
					if ($is_test_db) $trio_parameters .= " -db NGSD_TEST -user unknown";
					$parser->execTool("Tools/db_queue_analysis.php", $trio_parameters);
				}
			}
		}
		else trigger_error("Should not happen!", E_USER_ERROR);

	}
	else if ($queue_basecalling)
	{
		if ($force_basecalling) trigger_error("Full basecalling requested!", E_USER_NOTICE);
		else trigger_error("Basecall models do not match: On-device basecall model: '{$on_device_basecall_model}', requested basecall model: '{$basecall_model}'. Full (re-)basecalling is needed!", E_USER_WARNING);

		//move basecalled files
		$basecall_suffix = "";
		if (str_contains($on_device_basecall_model, "hac")) $basecall_suffix = "_hac";
		else if (str_contains($on_device_basecall_model, "sup")) $basecall_suffix = "_sup";
		else if (str_contains($on_device_basecall_model, "fast")) $basecall_suffix = "_fast";
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
		if (file_exists($out_bam)) trigger_error("Output BAM file '{$out_bam}' already exists, aborting.", E_USER_ERROR);
		

		//queue full re-basecalling
		$basecalling_cmd = create_basecall_command($run_dir, $out_bam, $basecall_model, $min_qscore, end($subdirs)."/bam_pass/", $queue_sample, false, $basecalling_copy_pod5);
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
		trigger_error("Basecall models do not match: On-device basecall model: '{$on_device_basecall_model}', requested basecall model: '{$basecall_model}'. Please check setting or re-do basecalling!", E_USER_ERROR);
	}


}

if ($queue_sample)
{
	// update sequencing run analysis status
	$hash = $db->prepare("UPDATE sequencing_run SET sequencing_run.status='analysis_started' WHERE sequencing_run.name = :run_name AND sequencing_run.status IN ('run_started', 'run_finished')");
	$db->bind($hash, "run_name", $run_name);
	$db->execute($hash, true);
}


?>