<?php 
/** 
	@page copy_ont_run
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$repo_folder = "/mnt/storage2/megSAP/pipeline"; //fixed absolute path to make the tests work for all users

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("copy_ont_run", "Copy ONT flow cell data (FASTQ, BAM).");
$parser->addString("run_name",  "Run name.", false);
$parser->addString("run_dir",  "Input data folder.", false);

$parser->addFlag("prefer_bam", "Prefer BAM over FASTQ, if available.");
$parser->addFlag("fastq", "Copy FASTQ files.");
$parser->addFlag("bam", "Copy BAM files (unaligned or aligned).");
$parser->addFlag("ignore_aligned", "Treat aligned BAM files as unaligned.");
$parser->addFlag("queue_sample", "Queue analysis of the sample.");

$parser->addString("single_fastq",  "Create single FASTQ file, without sample lookup and database checks.", true, "");

$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addInt("threads", "Number of threads to use for file merging and compression.", true, 4);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

$file_acccess_group = get_path("file_access_group", false);
if ($file_acccess_group == "") trigger_error("File access group not set in the settings.ini. File group will not be set!", E_USER_WARNING);

//absolute path
$run_dir = realpath($run_dir);
if (!file_exists($run_dir))
{
	trigger_error("Run directory '{$run_dir}' does not exists!", E_USER_ERROR);
}

//set ulimit
exec2("ulimit -n 10000");

//find subdirectory in run directory
$subdirs = array_values(array_diff(scandir($run_dir), array("..", ".")));
if (count($subdirs) === 0)
{
	trigger_error("No subdirectories found in '{$run_dir}'.", E_USER_ERROR);
}
foreach ($subdirs as &$subdir)
{
	$subdir = "{$run_dir}/{$subdir}";
}

//check for pod5_skip subdirectory (partially basecalled flow cell)
foreach ($subdirs as $subdir)
{
	if (is_dir("{$subdir}/pod5_skip"))
	{
		list($stdout, $stderr, $exit_code) = $parser->exec("du", "-b --summarize {$subdir}/pod5_skip");
		$folder_size = intval(explode("\t", $stdout[0])[0]);
		if ($folder_size > 5*1024*1024) // > 5MB
		{
			trigger_error("'pod5_skip' directory present in '{$subdir}', some data has not been basecalled!", E_USER_ERROR);
		}
		
	}
}

//check available data
function check_data_available($path, $subdir_name, $suffix)
{
	return file_exists($path) && is_dir($path) && (count(glob("{$path}/*/{$subdir_name}/*{$suffix}")) > 0);
}

$fastq_available = check_data_available($run_dir, "fastq_pass", ".fastq.gz") || check_data_available($run_dir, "fastq_pass/barcode??", ".fastq.gz");
if (!$fastq_available && ($fastq || $single_fastq))
{
	trigger_error("No FASTQ files available!", E_USER_ERROR);
}
$fastq_paths_glob = "{$run_dir}/*/fastq_pass";


if ($single_fastq !== "")
{
	trigger_error("Copy and merge FASTQ files to given output file.", E_USER_NOTICE);

	if (file_exists($single_fastq))
	{
		trigger_error("Output FASTQ file '{$single_fastq}' already exists, aborting.", E_USER_ERROR);
	};

	exec2("find {$fastq_paths_glob} -name '*.fastq.gz' -type f -exec zcat -f {} + | pigz -p {$threads} -c > {$single_fastq}");
	trigger_error("FASTQ saved in {$single_fastq}.", E_USER_NOTICE);
	exit();
}

$bam_available = check_data_available($run_dir, "bam_pass", ".bam") || check_data_available($run_dir, "bam_pass/barcode??", ".bam");
if (!$bam_available && $bam)
{
	trigger_error("No BAM files available!", E_USER_ERROR);
}

if (!$bam_available && !$fastq_available && $prefer_bam)
{
	trigger_error("No FASTQ/BAM files are available!", E_USER_ERROR);
}
$bam_paths_glob = "{$run_dir}/*/bam_pass";

//database connection
$db_con = DB::getInstance($db);

$result = $db_con->executeQuery("SELECT id, fcid FROM sequencing_run WHERE name = '{$run_name}'");
if (count($result) !== 1)
{
	trigger_error("Sequencing run with the name {$run_name} not found in database.", E_USER_ERROR);
}
$run_id = $result[0]['id'];
$flowcell_id = $result[0]['fcid'];

//check flowcell ID
foreach ($subdirs as $subdir)
{
	if (!mb_strpos(basename($subdir), $flowcell_id))
	{
		trigger_error("Flowcell ID '{$flowcell_id}' not found in directory name '{$run_dir}'.", E_USER_ERROR);
	}
}

//check random BAM file to find out which information is present
$bam_files = glob("{$bam_paths_glob}/*.bam");
if (count($bam_files) !== 0)
{
	$bam_random_file = $bam_files[array_rand($bam_files, 1)];

	$ret = $parser->exec(get_path("samtools"), "view --count --exclude-flags 0x900 {$bam_random_file}");
	$num_records = intval($ret[0][0]);
	$ret = $parser->exec(get_path("samtools"), "view --count --tag ML {$bam_random_file}");
	$num_base_mods = intval($ret[0][0]);
	$ret = $parser->exec(get_path("samtools"), "view --count --require-flags 0x4 {$bam_random_file}");
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

//import run QC
$parser->execTool("NGS/runqc_parser_ont.php", "-name '{$run_name}' -run_dir {$run_dir} -db {$db} -force");

//find sample entered in database
$result = $db_con->executeQuery(<<<SQL
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
		$bam_paths_glob = "{$run_dir}/*/bam_pass/{$barcode}";
		$fastq_paths_glob = "{$run_dir}/*/fastq_pass/{$barcode}";
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

	//perform GenLab import before analysis (skip on test run)
	if ($db!="NGSD_TEST")
	{
		trigger_error("Importing information from GenLab...", E_USER_NOTICE);
		$args = [];
		$args[] = "-ps {$sample}";
		$parser->exec(get_path("ngs-bits")."/NGSDImportGenlab", implode(" ", $args), true);
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
			$pipeline[] = [get_path("samtools"), "merge --reference {$genome} --threads {$threads} -b - -o {$out_bam}"];
			$parser->execPipeline($pipeline, "merge aligned BAM files");
			$parser->indexBam($out_bam, $threads);

		}
		else
		{
			$pipeline[] = [get_path("samtools"), "cat --threads {$threads} -o {$out_bam} -b -"]; //no reference required
			$parser->execPipeline($pipeline, "merge unaligned BAM files");
		}
	}

	if ($fastq || ($prefer_bam && !$bam_available))
	{
		trigger_error("Copy and merge FASTQ files.", E_USER_NOTICE);

		$out_fastq = "{$out_dir}{$sample}.fastq.gz";
		if (file_exists($out_fastq))
		{
			trigger_error("Output FASTQ file '{$out_fastq}' already exists, aborting.", E_USER_ERROR);
		};

		exec2("find {$fastq_paths_glob} -name '*.fastq.gz' -type f -exec cat  {} + > {$out_fastq}");
		trigger_error("FASTQ saved in {$out_fastq}.", E_USER_NOTICE);
	}

	//apply file access permissions
	$parser->exec("chmod", "-R 775 {$out_dir}");
	if ($file_acccess_group != "") $parser->exec("chgrp", "-R {$file_acccess_group} {$out_dir}");

	if ($queue_sample)
	{
		$parser->execTool("NGS/db_queue_analysis.php", "-samples {$sample} -type 'single sample' -db {$db}".(($db=="NGSD_TEST")?" -user unknown":""));
	}
}

if ($queue_sample)
{
	// update sequencing run analysis status
	$db_con->executeStmt("UPDATE sequencing_run SET sequencing_run.status='analysis_started' WHERE sequencing_run.name = '{$run_name}' AND sequencing_run.status IN ('run_started', 'run_finished')");
}

?>