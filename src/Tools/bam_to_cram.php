<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("bam_to_cram", "Convert BAM to CRAM file.");
$parser->addInfile("bam", "Input BAM file.", false);
$parser->addOutFile("cram", "Output CRAM file.", false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addInt("threads", "Number of threads to use.", true, 1);
$parser->addEnum("format", "Output CRAM format", true, ["3.0", "3.1"], "3.0");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//convert BAM to CRAM
$parser->execApptainer("samtools", "samtools view", "-@ {$threads} -C --output-fmt-option version={$format} -T {$genome} -o {$cram} {$bam}", [$genome, $bam], [dirname($cram)]);

//index CRAM
$parser->execApptainer("samtools", "samtools index", "-@ {$threads} {$cram}", [dirname($cram)]);

?>