<?php
/** 
	@page delete_fastq
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);


//parse command line arguments
$parser = new ToolBase("delete_fastq", "Deletes all FASTQ files of a given processing system.");
$parser->addString("processing_system", "Processing system of the samples.", false);
$parser->addFlag("delete", "Deletes the FASTQ files. Warning: No additional warning before deleting files!");
extract($parser->parse($argv));

$ngsbits = get_path("ngs-bits");

//get samples with matching processing system
$sample_table_file = $parser->tempFile(".tsv");
$pipeline= [
	["{$ngsbits}NGSDExportSamples", "-system $processing_system -add_path SAMPLE_FOLDER"],
	["{$ngsbits}TsvFilter", "-filter 'system_name_short is $processing_system'"], //this is necessary because NGSDExportSamples performs fuzzy match
	["{$ngsbits}TsvSlice", "-cols name,path -out $sample_table_file"],
];
$parser->execPipeline($pipeline, "NGSD sample extraction", true);
$sample_table = load_tsv($sample_table_file);

//process samples
print "#ps_name\tfolder\tfastq_file_size\tbam_file_size\taction\terrors\n";
foreach($sample_table as $sample)
{	
	//init
	$ps_name = $sample[0];
	$s_name = substr($sample[0], 0, -3);
	$folder = $sample[1];
	$size = "0.0";
	$errors = array();
	
	// get FASTQ files:
	$files1 = glob($folder."/".$s_name."*_R1_00?.fastq.gz");
	$files2 = glob($folder."/".$s_name."*_R2_00?.fastq.gz");
	$fastq_files = array_merge($files1, $files2);
	
	// check FASTQ files
	if(count($files1)!=count($files2))
	{
		$errors[] = "Mismatching forward and reverse FASTQ count";
	}
	if(count($fastq_files)==0)
	{
		$errors[] = "No FASTQ files found";
	}

	// determine FASTQ file(s) size
	$fastq_file_size = 0;
	foreach($fastq_files as $fq_file)
	{
		$fastq_file_size += filesize($fq_file);
	}
	$size = number_format($fastq_file_size/1024/1024/1024, 3, ".", "");


	//check BAM file
	$bam_file = $folder."/".$ps_name.".bam";
	$bam_exists = file_exists($bam_file) && file_exists($bam_file.".bai"); 
	if(!$bam_exists)
	{
		$errors[] = "No BAM/BAI file found!";
	}
	
	//check BAM file size
	$bam_file_size = 0;
	if ($bam_exists && count($fastq_files)>0)
	{
		$bam_file_size = filesize($bam_file);
		
		if ($bam_file_size / $fastq_file_size < 0.4)
		{
			$errors[] = "BAM file smaller than 40% of FASTQ";
		}
	}
	
	//determine action
	$action = "SKIP";
	if (count($errors)==0)
	{
		$action = "DELETE ".count($fastq_files)." files".($delete ? "" : " (dry-run)");
	}
	
	//output
	print "{$ps_name}\t{$folder}\t".number_format($fastq_file_size/1024/1024/1024, 3, ".", "")."\t".number_format($bam_file_size/1024/1024/1024, 3, ".", "")."\t{$action}\t".implode(" ", $errors)."\n";
	
	//delete
	if ($delete && starts_with($action, "DELETE"))
	{
		foreach($fastq_files as $fq_file)
		{
			unlink($fq_file);
		}
	}
}

?>