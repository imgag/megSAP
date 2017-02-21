<?php

/*
	@page indel_realign
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("indel_realign", "Perform indel realignment using GATK.");
$parser->addInfile("in",  "Input BAM file.", false);
$parser->addOutfile("out", "Output BAM file.", false);
//optional
$parser->addString("genome", "Path to GATK reference genome FASTA file.", true, get_path("data_folder")."genomes/GATK/hg19.fa");
extract($parser->parse($argv));

//execute indel realignment step
$intervalsFile = $parser->tempFile(".intervals");
$parser->exec(get_path("GATK"), "-T RealignerTargetCreator -I $in -R $genome -o $intervalsFile", true);
$parser->exec(get_path("GATK"), "-T IndelRealigner -I $in -R $genome --targetIntervals $intervalsFile -o $out -rf NotPrimaryAlignment", true);

//rename .BAI file to .BAM.BAI
rename(substr($out, 0, -4).".bai", $out.".bai"); 
?>