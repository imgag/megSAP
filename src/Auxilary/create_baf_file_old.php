<?php
/** 
	@page create_baf_file_old
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("create_baf_file_old", "Creates BAF file from GSVar file and bam file.");
$parser->addOutfile("out_file", "Target to store BAF file.", false);
$parser->addInfile("gsvar", "Input GSvar file.", false);
$parser->addInfile("bam", "Input BAM file.", false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));


//Abort if out_file exists to prevent interference with other jobs
if(file_exists($out_file)) return;
$genome = genome_fasta($build);

$tmp_out = $parser->tempFile(".tsv");
$parser->execApptainer("ngs-bits", "VariantAnnotateFrequency", "-in {$gsvar} -bam {$bam} -depth -out {$tmp_out} -ref {$genome}", [$gsvar, $bam, $genome]);

$in_handle  = fopen2($tmp_out,"r");
$out_handle = fopen2($out_file,"w");

while(!feof($in_handle))
{

    $line = trim(fgets($in_handle));
    if(starts_with($line,"#")) continue;
    if(empty($line)) continue;
    
    $parts = explode("\t",$line);
    list($chr,$start,$end,$ref,$alt) = $parts;
    if(strlen($ref) != 1 || strlen($alt) != 1 ) continue; //Skip INDELs

    $freq = $parts[count($parts)-2]; //freq annotation ist 2nd last col
    $depth= $parts[count($parts)-1]; //depth annotation is last col
    fputs($out_handle,"{$chr}\t{$start}\t{$end}\t{$chr}_{$start}\t{$freq}\t{$depth}\n");
}
fclose($in_handle);
fclose($out_handle);

?>