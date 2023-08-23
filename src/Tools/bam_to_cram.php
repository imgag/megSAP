<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("bam_to_cram", "Convert sample BAM to CRAM for GRCh38.");
$parser->addString("ps", "Processed sample name.", false);
$parser->addInt("threads", "Number of threads to use.", true, 1);
extract($parser->parse($argv));

//init
$db = DB::getInstance("NGSD");

//get BAM
$info = get_processed_sample_info($db, $ps);
$ps_folder = $info['ps_folder'];
$bam = $ps_folder."/".$ps.".bam";

if (!file_exists($bam))
{
	print "$ps skipped - no BAM found!\n";
}
else
{
	//check that BAM is not mapped on genome with false duplication maked (the region is in the center of the biggest false duplication, so no read should be mapped there if the genome is correctly masked)
	$samtools = get_path("samtools");
	list($stdout) = exec2("{$samtools} view {$bam} chr21:5966593-6161371 | wc -l");
	$mapped_read_count = trim(implode("", $stdout));
	if ($mapped_read_count>0)
	{
		print "$ps skipped - BAM was produced with GRCh38 where false duplication are not masked!\n";
	}
	else
	{
		//convert BAM to CRAM
		$genome = genome_fasta("GRCh38");
		$cram = $ps_folder."/".$ps.".cram";
		exec2("{$samtools} view -@ {$threads} -C -T {$genome} -o {$cram} {$bam}");
		
		//index CRAM
		exec2("{$samtools} index -@ {$threads} {$cram}");
		
		//delete BAM and BAI
		exec2("rm -f {$bam}");
		exec2("rm -f {$bam}.bai");
	}
}

?>