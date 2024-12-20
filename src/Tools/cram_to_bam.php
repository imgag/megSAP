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
$genome = genome_fasta($build);

//convert CRAM to BAM
$parser->execApptainer("samtools", "samtools view", "-@ {$threads} -b -T {$genome} -o {$bam} {$cram}", [$genome, $cram], [dirname($bam)]);

//index BAM
$parser->execApptainer("samtools", "samtools index", "-@ {$threads} {$bam}", [], [dirname($bam)]);

?>