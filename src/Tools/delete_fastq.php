<?php
/** 
	@page delete_fastq
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);


//parse command line arguments
$parser = new ToolBase("delete_fastq", "Deletes all FASTQ files of a given processing system.");
$parser->addString("processing_system", "Processing system of the FASTQ files.", false);
$parser->addFlag("delete", "Deletes the FASTQ files. Warning: No additionaol warning before deleting files!");
extract($parser->parse($argv));

$ngsbits = get_path("ngs-bits");


// get all samples with matching processing system
$sample_table_file = $parser->tempFile(".tsv");
$parser->exec("{$ngsbits}NGSDExportSamples", "-out $sample_table_file -system $processing_system -add_path", true);
$sample_table = load_tsv($sample_table_file);

// overall stats:
$overall_fastq_size = 0;


foreach($sample_table as $sample)
{
	// print sample infos
	print substr($sample[0], 0, -3)." at: ".$sample[18]."\n";

	// get FASTQ files:
	$in_for = $sample[18].substr($sample[0], 0, -3)."*_R1_00?.fastq.gz";
	$in_rev = $sample[18].substr($sample[0], 0, -3)."*_R2_00?.fastq.gz";
	$files1 = glob($in_for);
	$files2 = glob($in_rev);

	// check if forward and reverse file count is equal
	if(count($files1)!=count($files2))
	{
		trigger_error("Found mismatching forward and reverse read file count!\n Forward: $in_for\n Reverse: $in_rev.", E_USER_WARNING);
		print "\033[0;31m Mismatching forward and reverse read file count! Skipping deletion. \033[0m\n\n\n";
		continue;
	}
	if(count($files1) == 0)
	{
		trigger_error("No FASTQ files found!", E_USER_WARNING);
		
		print "\033[0;31m No FASTQ files! Skipping deletion. \033[0m\n\n\n";
		continue;
	}

	// print file names and calculate FASTQ file size:
	$fastq_file_size = 0;
	print "FASTQ files:\n";
	foreach($files1 as $fq_file)
	{
		$fastq_file_size += filesize($fq_file);
		print "\t".$fq_file."\n";
	}
	foreach($files2 as $fq_file)
	{
		$fastq_file_size += filesize($fq_file);
		print "\t".$fq_file."\n";
	}
	
	// print FASTQ filesize
	print "FASTQ file size:\t\t ".number_format($fastq_file_size/1024/1024/1024, 3)." GB \n";


	// check if BAM file exists
	$bam_file = $sample[18].$sample[0].".bam";
	if(!file_exists($bam_file) || !file_exists($bam_file.".bai"))
	{
		trigger_error("No BAM file found!", E_USER_WARNING);
		print "\033[0;31m No BAM file found! Skipping deletion. \033[0m\n\n\n";
		continue;
	}

	// print BAM file and file size
	$bam_file_size = filesize($bam_file);
	print "\nBAM file:\n\t".$bam_file."\n";
	print "BAM file size:\t\t\t ".number_format($bam_file_size/1024/1024/1024, 3)." GB \n";


	// check if BAM file size fits 
	$bam_fastq_ratio = $bam_file_size / $fastq_file_size;

	if ($bam_fastq_ratio > 0.5)
	{
		print "\033[0;32m Bam file larger than 50% of FASTQ: FASTQ files ready for deletion. \033[0m\n";

		$overall_fastq_size += $fastq_file_size;

		if ($delete)
		{
			print "\tDeleting files...";
			foreach($files1 as $fq_file)
			{
				unlink($fq_file);
			}
			foreach($files2 as $fq_file)
			{
				unlink($fq_file);
			}
		}
		else
		{
			print "\tSimulation run. No files will be deleted.";
		}
	}
	else
	{
		print "\033[0;31m Bam file smaller than 50% of FASTQ: FASTQ will not be deleted. \033[0m\n\n\n";
		continue;
	}


	print "\n\n\n";
}

print "\nOverall size of deleted files:\t\t".number_format($overall_fastq_size/1024/1024/1024, 0)." GB\n";

?>