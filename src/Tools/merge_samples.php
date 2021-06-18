<?php
/** 
	@page merge_samples 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("merge_samples", "Merges a processed sample into another processed sample.");
$parser->addString("ps", "Processed sample to merged into the second sample.", false);
$parser->addString("into", "Processed sample into which the first sample is merged.", false);
extract($parser->parse($argv));

//get processed sample infos
$db = DB::getInstance("NGSD");
$info1 = get_processed_sample_info($db,$ps);
$info2 = get_processed_sample_info($db,$into);

$ngsbits = get_path("ngs-bits");

//check input folders
$folder1 = $info1['ps_folder'];
if (!is_dir($folder1))
{
	trigger_error("Sample folder '$folder1' does not exist!", E_USER_ERROR);
}
$folder2 = $info2['ps_folder'];
if (!is_dir($folder2))
{
	trigger_error("Sample folder '$folder2' does not exist!", E_USER_ERROR);
}

//make sure the "+merged_samples" folder exists
$backup_folder = dirname($folder1)."/+merged_samples";
if (!file_exists($backup_folder))
{
	if (!mkdir($backup_folder))
	{
		trigger_error("Could not create backup folder '$backup_folder'!", E_USER_ERROR);
	}
}

//check if FASTQ files or BAM in target folder exist
$target_fastq_files = glob($folder2."/*.fastq.gz");
$target_bam_file = $folder2."/".$into.".bam";
if (count($target_fastq_files) > 0)
{
	//FASTQ files found -> nothing to do
}
elseif (file_exists($target_bam_file)) 
{
	trigger_error("No FASTQ files found in Sample folder of '${into}'. Generating FASTQs from BAM.", E_USER_WARNING);
	// recreate FASTQ files from BAM
	$in_fq_for = $folder2."/${into}_BamToFastq_R1_001.fastq.gz";
	$in_fq_rev = $folder2."/${into}_BamToFastq_R2_001.fastq.gz";
	$tmp1 = $parser->tempFile(".fastq.gz");
	$tmp2 = $parser->tempFile(".fastq.gz");
	$parser->exec("{$ngsbits}BamToFastq", "-in ${target_bam_file} -out1 $tmp1 -out2 $tmp2", true);
	$parser->moveFile($tmp1, $in_fq_for);
	$parser->moveFile($tmp2, $in_fq_rev);
}
else
{
	trigger_error("Could not find any sequencing data (neither FASTQ nor BAM) in Sample folder of '$into'! Abort merging.", E_USER_ERROR);
}

//check if FASTQ files or BAM in source folder exist
$source_fastq_files = glob($folder1."/*.fastq.gz");
$source_bam_file = $folder1."/".$ps.".bam";
if (count($source_fastq_files) > 0)
{
	//FASTQ files found -> move them to the target folder
	exec2("mv $folder1/*.fastq.gz $folder2/");
}
elseif (file_exists($source_bam_file)) 
{
	trigger_error("No FASTQ files found in Sample folder of '${ps}'. Generating FASTQs from BAM.", E_USER_WARNING);
	// recreate FASTQ files from BAM
	$in_fq_for = $folder2."/${ps}_BamToFastq_R1_001.fastq.gz";
	$in_fq_rev = $folder2."/${ps}_BamToFastq_R2_001.fastq.gz";
	$tmp1 = $parser->tempFile(".fastq.gz");
	$tmp2 = $parser->tempFile(".fastq.gz");
	$parser->exec("{$ngsbits}BamToFastq", "-in ${source_bam_file} -out1 $tmp1 -out2 $tmp2", true);
	$parser->moveFile($tmp1, $in_fq_for);
	$parser->moveFile($tmp2, $in_fq_rev);

	//delete BAM
	unlink($source_bam_file);	
}
else
{
	trigger_error("Could not find any sequencing data (neither FASTQ nor BAM) in Sample folder of '$ps'! Abort merging.", E_USER_ERROR);
}

// backup all other files in target folder
exec2("mv $folder1 $backup_folder/");


//NGSD: remove variants/qc for 'ps'
$ps_id = $info1['ps_id'];
$db->executeStmt("DELETE FROM detected_variant WHERE processed_sample_id='$ps_id'");
$db->executeStmt("DELETE FROM processed_sample_qc WHERE processed_sample_id='$ps_id'");

//NGSD: mark samples as merged
$db->executeStmt("INSERT INTO `merged_processed_samples` (`processed_sample_id`, `merged_into`) VALUES ({$ps_id},".$info2['ps_id'].")");

?>