<?php
/**
	@page detect_msi
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$parser = new ToolBase("detect_msi", "Detects microsatellite instabilities via MSIsensor-pro");
$parser->addInfile("n_bam","BAM file containing normal data",true);
$parser->addInfile("t_bam","BAM file containing tumor data",true);
$parser->addInfile("msi_ref","MSIsensor-pro scan file for the genome used in the samples.",true);
$parser->addString("out","Path to the output file",true);
//optional
$parser->addString("build", "The reference genome build to use. ", true, "GRCh38");
$parser->addFlag("keep_status_files","Keep MSI status files");
$parser->addInt("threads", "The maximum number of threads to use.", true, 1);

extract($parser->parse($argv));

$parameters = 	" msi -b $threads -d $msi_ref -n $n_bam -t $t_bam -o $out ";
$parser->exec(get_path("msisensor"), $parameters);

if (! $keep_status_files)
{
	if (file_exists($out."_dis")) exec2("rm {$out}_dis");
	if (file_exists($out."_germline")) exec2("rm {$out}_germline");
	if (file_exists($out."_somatic")) exec2("rm {$out}_somatic");
}

if (! file_exists($out))
{
	trigger_error("MSI calling did not exit successfully. No MSIsensor output file was created.",E_USER_ERROR);
}

//prepend # to the first line:
file_put_contents($out, "#".file_get_contents($out));
?>