<?php

/**
	@page tool_stub
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_control_freec", "CNV detection with Control-FREEC.");

// mandatory arguments
$parser->addInfile("tumor_bam", "Tumor BAM file.", false);
$parser->addInfile("normal_bam", "Normal BAM file.", false);
$parser->addString("out_dir", "Output directory.", false);

// optional arguments

// extract arguments
extract($parser->parse($argv));

if (!is_dir($out_dir))
{
	mkdir($out_dir);
}

// Control-FREEC configuration
$config = <<<CONFIG
[general]
chrLenFile = /mnt/users/ahmattj1/software/upstream/FREEC/GRCh37.fa.fai
window = 0
ploidy = 2,3
outputDir = {$out_dir}

#sex=XY
breakPointType=4
chrFiles = /mnt/users/ahmattj1/software/upstream/FREEC/chr-files

maxThreads=10

breakPointThreshold=1.2
noisyData=TRUE
printNA=FALSE

readCountThreshold=50

#consider this
#contaminationAdjustment=TRUE

degree=1

[sample]
mateFile = {$tumor_bam}
inputFormat = BAM
mateOrientation = FR

[control]
mateFile = {$normal_bam}
inputFormat = BAM
mateOrientation = FR

[BAF]

[target]
captureRegions = /mnt/share/data/enrichment/ssHAEv6_2017_01_05.bed
CONFIG;


$conf = $parser->tempFile();
file_put_contents($conf, $config);
$parser->log("Config file:", file($conf));

//run command
$parser->exec("/mnt/users/ahmattj1/software/upstream/FREEC/src/freec",
	"-conf {$conf}", true);