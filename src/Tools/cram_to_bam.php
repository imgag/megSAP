<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("cram_to_bam", "Convert CRAM to BAM file.");
$parser->addInfile("cram", "Input CRAM file.", false);
$parser->addOutFile("bam", "Output BAM file.", false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addInt("threads", "Number of threads to use.", true, 1);
extract($parser->parse($argv));

//init
$samtools = get_path("samtools");
$genome = genome_fasta($build);

//check that CRAM is mapped on genome with false duplication masked (the region is in the center of the biggest false duplication, so no read should be mapped there if the genome is correctly masked)
if ($build=="GRCh38")
{
	list($stdout) = exec2("{$samtools} view -T {$genome} {$cram} chr21:5966593-6161371 | wc -l");
	$mapped_read_count = trim(implode("", $stdout));
	if ($mapped_read_count>0)
	{
		trigger_error("CRAM file '$cram' cannot be converted to BAM - it was produced with GRCh38 where false duplication are not masked!", E_USER_ERROR);
	}
}

//convert CRAM to BAM
exec2("{$samtools} view -@ {$threads} -b -T {$genome} -o {$bam} {$cram}");

//index BAM
exec2("{$samtools} index -@ {$threads} {$bam}");

?>