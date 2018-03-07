<?php
/**
	@page detect_msi
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("detect_msi", "Detects microsatellite instabilities via MANTIS");
$parser->addInfile("n_bam","Bam file containing normal data",true);
$parser->addInfile("t_bam","Bam file containing tumor data",true);
$parser->addInt("threads","Maximum number of CPU threads to use",true,1);
$parser->addOutfile("out","Path to the output file",true);
$parser->addInfile("bed_file","Bed file that contains target region",true);
//optional
$parser->addString("build", "The reference genome build to use. ", true, "GRCh37");
extract($parser->parse($argv));

$parameters = "-n $n_bam -t $t_bam  -b $bed_file -o $out --threads $threads --genome $build";

$parser->exec(get_path("mantis"),$parameters,true,true);




?>