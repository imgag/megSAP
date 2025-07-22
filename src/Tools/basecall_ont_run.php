<?php 
/** 
	@page basecall_ont_run
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("basecall_ont_run", "Performs dorado basecalling on ONT run folder.");
$parser->addInfile("run_dir",  "Input data folder.", false);
$parser->addOutfile("out_bam", "Output file path to store (basecalled) unmapped BAM.", false);

//optional
$parser->addString("basecall_model", "Model used for base calling: hac=high accuracy, sup=super accuracy", true, "hac");
$parser->addInt("min_qscore", "Minimal read qScore for basecalling", true, 9);
$parser->addFlag("queue_sample", "Queue analysis of the sample.");
$parser->addFlag("skipped_only", "Use already basecalled BAMs and recall only skipped POD5s.");
$parser->addFlag("copy_pod5", "Copy pod5 files to /tmp to reduce IO/network load.");
$parser->addFlag("file_based_calling", "Call each POD5 file separately.");

extract($parser->parse($argv));


//check GPU server


//init
//check basecall model
if ($basecall_model == "hac")
{
	$basecall_model = "hac@".get_path("dorado_model_version").",5mCG_5hmCG";
}
else if ($basecall_model == "sup")
{
	$basecall_model = "sup@".get_path("dorado_model_version").",5mCG_5hmCG";
}
else
{
	trigger_error("Custom basecall model '{$basecall_model}' provided!", E_USER_WARNING);
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


$bam_files = array();
if ($skipped_only)
{
	//get BAMs
	$bam_files = glob("{$run_dir}/*/bam_pass/*.bam");

	if (count($bam_files) < 1) trigger_error("No BAM files found in run folder!", E_USER_ERROR);

	//check model used for basecalling of bam files (test 10 random files)
	for ($i=0; $i < 10; $i++) 
	{ 
		$idx = random_int(0, count($bam_files));
		$model = get_basecall_model($bam_files[$idx]);

		if (starts_with($basecall_model, "hac") && !starts_with($model, "hac")) trigger_error("Run should be basecalled with high accuracy, but already basecalled BAMs are not!" , E_USER_ERROR);
		if (starts_with($basecall_model, "sup") && !starts_with($model, "sup")) trigger_error("Run should be basecalled with super accuracy, but already basecalled BAMs are not!" , E_USER_ERROR);
	}

}


//get POD5 location
$pod5_location_pattern = "{$run_dir}/*/pod5*/";
if ($skipped_only) $pod5_location_pattern = "{$run_dir}/*/pod5_skip/";

//copy POD5s
$pod5_files = glob($pod5_location_pattern."*.pod5");
if ($copy_pod5)
{
	$pod5_temp = $parser->tempFolder("pod5_temp");
	foreach($pod5_files as $file)
	{
		$parser->copyFile($file, $pod5_temp."/".basename($file));
	}

	$pod5_location = [$pod5_temp];
	$pod5_files = glob($pod5_temp."/*.pod5"); //use copied pod5s 
}
else
{
	$pod5_location = glob($pod5_location_pattern);
}

if (count($pod5_location) < 1) trigger_error("No POD5 files found in '{$pod5_location_pattern}'!", E_USER_ERROR);


//TODO: create container
// perform basecalling
$bams_to_merge = [];
$dorado_model_path = get_path("dorado_model_path");

if ($file_based_calling)
{
	foreach ($pod5_files as $pod5_file) 
	{
		$tmp_bam = $parser->tempFile(basename2($pod5_file).".mod.unmapped.bam");
		$parser->exec(get_path("dorado"), "basecaller --models-directory {$dorado_model_path} {$basecall_model} {$pod5_file} --min-qscore {$min_qscore} > {$tmp_bam}");
		$bams_to_merge[] = $tmp_bam;
	}
}
else
{
	foreach ($pod5_location as $folder) 
	{
		$tmp_bam = $parser->tempFile(".mod.unmapped.bam");
		$parser->exec(get_path("dorado"), "basecaller --models-directory {$dorado_model_path} -r {$basecall_model} {$folder} --min-qscore {$min_qscore} > {$tmp_bam}");
		$bams_to_merge[] = $tmp_bam;
	}
}




if ($skipped_only)
{
	//add already basecalled bams
	list($stdout, $stderr, $ec) = $parser->exec("find", "{$run_dir}/*/bam_pass/*.bam -name '*.bam' -type f");

	$bams_to_merge = array_merge($bams_to_merge, $stdout);
}


if (count($bams_to_merge) > 1)
{
	//merge bams and copy to output location
	$merged_bam = $parser->tempFile(".merged.mod.unmapped.bam");
	$bam_list = $parser->tempFile(".bams_to_merge.txt");
	file_put_contents($bam_list, implode("\n", $bams_to_merge));
	$dorado_model_path = get_path("dorado_model_path");
	$parser->execApptainer("samtools", "samtools cat", "--threads 20 -o {$merged_bam} -b {$bam_list}", [$run_dir], []);

	
	//copy output
	$parser->copyFile($merged_bam, $out_bam);
}
else
{
	//only copy necessary
	$parser->copyFile($bams_to_merge[0], $out_bam);
}



//queue sample
if ($queue_sample)
{
	$sample_name = basename2($out_bam);
	$parser->execTool("Tools/db_queue_analysis.php", "-samples {$sample_name} -type 'single sample'");
}



?>