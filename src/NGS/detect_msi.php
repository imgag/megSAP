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
$parser->addString("build", "The reference genome build to use. ", true, "GRCh38");
$parser->addFlag("keep_status_files","Keep MSI status files");
$parser->addFlag("is_exome","Use standard parameters for MSI calling dedicated for exomes.");
extract($parser->parse($argv));

$output_folder = realpath(dirname($out));
$output_path = $output_folder . '/' . 'mantis';
$parameters = 	"-n $n_bam -t $t_bam  -b $bed_file -o $output_path --threads $threads --genome ".genome_fasta($build);
if($is_exome) //use non-standard MANTIS-parameters for whole exomes recommended by MANTIS
{
	$parameters .=  " -mrq 20 -mlq 25 -mlc 20 -mrr 1";
}

$parser->exec(get_path("python27"), get_path("mantis")." ".$parameters,true,false);

//adds comment char "#" to certain lines of the MANTIS status file
function parse_mantis_status_file($out_file_name)
{	
	$file = fopen2($out_file_name,'r');
	$lines = array();
	$i = 0;
	while(!feof($file))
	{
		if($i == 0 || $i >= 4) //add comment char in line 0 and in lines 4-7
		{
			$lines[] = '#' . fgets($file);
		}
		else
		{
			$lines[] = fgets($file);
		}
		++$i;
	}
	fclose($file);
	
	$content_as_string = implode($lines);
	file_put_contents($out_file_name,$content_as_string);
}

//main

//remove not needed MANTIS output files
if(!$keep_status_files)
{
	if(file_exists($output_folder ."/mantis.kmer_counts")) exec2("rm $output_folder"."/mantis.kmer_counts");
	if(file_exists($output_folder ."/mantis.kmer_counts_filtered")) exec2("rm $output_folder"."/mantis.kmer_counts_filtered");
	if(file_exists($output_folder . "/mantis")) exec2("rm $output_path");
}

//rename important main output file
if(file_exists($output_folder ."/mantis.status"))
{
	exec2("mv $output_path".".status"." $out");
}
else
{
	trigger_error("MSI calling did not exit successfully. No MANTIS output file was created.",E_USER_WARNING);
}

parse_mantis_status_file($out);
?>