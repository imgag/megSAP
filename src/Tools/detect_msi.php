<?php
/**
	@page detect_msi
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$parser = new ToolBase("detect_msi", "Detects microsatellite instabilities via MSIsensor-pro");
$parser->addInfile("n_bam", "BAM/CRAM file containing normal data",false);
$parser->addInfile("t_bam", "BAM/CRAM file containing tumor data",false);
$parser->addString("msi_ref", "MSIsensor-pro scan file for the genome used in the samples.",false);
$parser->addString("out", "Path to the output file",false);
//optional
$parser->addInfile("ref", "reference genome. Required for CRAM input, but will be autofilled if possible using the build.", true);
$parser->addString("build", "The reference genome build to use. ", true, "GRCh38");
$parser->addFlag("keep_status_files", "Keep MSI status files");
$parser->addInt("threads", "The maximum number of threads to use.", true, 1);

extract($parser->parse($argv));

if ($ref == "")
{
	$ref = genome_fasta($build);
}

if(!file_exists($msi_ref)) // create msi-ref file
{
	print("Could not find loci reference file $msi_ref. Trying to generate it.\n");
	$parser->execApptainer("msisensor-pro", "msisensor-pro-v1.3.0", "scan -d $ref -o $msi_ref", [$ref, $msi_ref], [dirname($ref), dirname($msi_ref)]);

	//remove sites with more than 100 repeat_times as that crashes dragen MSI:
	$out_lines = [];
	foreach(file($msi_ref) as $line)
	{
		#$chr, $pos, $repeat_unit_len, $repeat_unit_bin, $repeat_times, ...
		$parts = explode("\t", $line);

		if ($parts[0] != "chromosome" && !chr_check($parts[0], 22, false)) continue;
		if (is_numeric($parts[4]) && floatval($parts[4]) > 100) continue;
		$out_lines[] = $line;
	}

	file_put_contents($msi_ref, implode("", $out_lines));
}

$parameters = 	"msi -b $threads -d $msi_ref -n $n_bam -g $ref -t $t_bam -o $out";

$parser->log("MSI ref file: $msi_ref");

$parser->execApptainer("msisensor-pro", "msisensor-pro-v1.3.0", $parameters, [$msi_ref, $n_bam, $t_bam], [dirname($out)]);

if (! $keep_status_files)
{
	if (file_exists($out."_dis")) exec2("rm {$out}_dis");
	if (file_exists($out."_all")) exec2("rm {$out}_all");
	if (file_exists($out."_unstable")) exec2("rm {$out}_unstable");
}

if (! file_exists($out))
{
	trigger_error("MSI calling did not exit successfully. No MSIsensor output file was created.",E_USER_ERROR);
}

//prepend # to the first line:
file_put_contents($out, "#".file_get_contents($out));
?>