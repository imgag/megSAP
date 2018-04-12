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
$parser->addFlag("keep_status_files","Keep MSI status files");
extract($parser->parse($argv));

$output_folder = realpath(dirname($out));
$output_path = $output_folder . '/' . 'mantis';
$parameters = "-n $n_bam -t $t_bam  -b $bed_file -o $output_path --threads $threads --genome $build";
$parser->exec(get_path("mantis"),$parameters,true,false);

//remove status files
if(!$keep_status_files)
{
	if(file_exists($output_folder ."/mantis.kmer_counts")) exec2("rm $output_folder"."/mantis.kmer_counts");
	if(file_exists($output_folder ."/mantis.kmer_counts_filtered")) exec2("rm $output_folder"."/mantis.kmer_counts_filtered");
	if(file_exists($output_folder ."/mantis.status")) exec2("rm $output_folder"."/mantis.status");
}

//rename main output file
if(file_exists($output_folder . "/mantis")) exec2("mv $output_path $out");


?>