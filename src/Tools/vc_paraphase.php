<?php 
/** 
	@page vc_paraphase
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_paraphase", "Caller for highly similar paralogous genes.");
$parser->addString("folder", "Sample folder for output files.", false);
$parser->addString("name", "Base file name, typically the processed name ID (e.g. 'GS120001_01').", false);

//optional
$parser->addInfile("local_bam",  "Optional alternative local BAM file (otherwise BAM from folder will be used).", false);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");

extract($parser->parse($argv));

//init
$genome = genome_fasta($build);
if ($local_bam != "")
{
	// use provided BAM
	$bam = $local_bam;
}
else
{
	//extract from folder
	if (file_exists("{$folder}/{$name}.bam")) $bam = "{$folder}/{$name}.bam";
	else if (file_exists("{$folder}/{$name}.cram")) $bam = "{$folder}/{$name}.cram";
	else trigger_error("No BAM/CRAM file found in output folder!", E_USER_ERROR);
}

//create output folder
$output_folder = "{$folder}/paraphase"; 
if (!file_exists($output_folder)) mkdir($output_folder);

//run paraphase
$args = [];
$args[] = "-p {$name}";
$args[] = "-b {$bam}";
$args[] = "-o {$output_folder}";
$args[] = "-r {$genome}";
$args[] = "-t {$threads}";

//set bind paths for modkit
$in_files = array($bam, $genome);
$out_files = array($output_folder);


$parser->execApptainer("paraphase", "paraphase", implode(" ", $args), $in_files, $out_files);

//TODO: Postprocessing


?>
