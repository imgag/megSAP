<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("bam_to_cram", "Convert BAM to CRAM file.");
$parser->addInfile("bam", "Input BAM file.", false);
$parser->addOutFile("cram", "Output CRAM file.", false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addInt("threads", "Number of threads to use.", true, 1);
extract($parser->parse($argv));

//init
$samtools = get_path("samtools");
$genome = genome_fasta($build);

//convert BAM to CRAM
exec2("{$samtools} view -@ {$threads} -C -T {$genome} -o {$cram} {$bam}");

//index CRAM
exec2("{$samtools} index -@ {$threads} {$cram}");

?>