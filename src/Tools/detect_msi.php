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
$parser->addInfile("target", "filters the msi_ref file to only include sites covered by the target region.", true);
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

if ($target != "")
{
	$parser->log("started filtering msi ref file.");
	$tmp_targeted_msi_ref = $parser->tempFile("detect_msi_targeted_msi_ref");
	$tmp_msi_ref_bed = $parser->tempFile("detect_msi_msi_ref_bed");
	$tmp_msi_ref_bed_filtered = $parser->tempFile("detect_msi_msi_ref_bed_filtered");
	
	//create tmp bed file from site start positions:
	$bed_lines = [];
	foreach(file($msi_ref) as $line)
	{
		$parts = explode("\t", $line);
		if($parts[0] == "chromosome") continue;
		$bed_lines[] = implode("\t", [$parts[0], $parts[1], intval($parts[1])+1]);
	}
	file_put_contents($tmp_msi_ref_bed, implode("\n", $bed_lines));
	file_put_contents("/mnt/storage2/users/ahott1a1/msi_ref_bed_test.bed", implode("\n", $bed_lines));
	
	//filter tmp bed by target:
	$parser->execApptainer("ngs-bits", "BedIntersect", "-in {$tmp_msi_ref_bed} -in2 {$target} -out {$tmp_msi_ref_bed_filtered}", [$tmp_msi_ref_bed, $target], [dirname($tmp_msi_ref_bed_filtered)]);
	exec("cp $tmp_msi_ref_bed_filtered /mnt/storage2/users/ahott1a1/msi_ref_bed_filtered_test.bed");
	$parser->log("msi ref exists? {$msi_ref} - ".file_exists($msi_ref));
	$handle_msi_ref = fopen2($msi_ref, "r");
	$handle_msi_ref_bed_filtered = fopen2($tmp_msi_ref_bed_filtered, "r");
	
	$sites_filtered = [];
	
	$sites_filtered[] = fgets($handle_msi_ref); //keep header
	$bed_filtered_line = fgets($handle_msi_ref_bed_filtered);
	if ($bed_filtered_line !== false)
	{
		list($bed_chr, $bed_start, $bed_end) = explode("\t", $bed_filtered_line);
	
		while(! feof($handle_msi_ref))
		{
			$msi_ref_line = fgets($handle_msi_ref);
			if (nl_trim($msi_ref_line) == "") continue; //skip empty lines
			
			$parts_msi_ref_line = explode("\t", $msi_ref_line);
			
			if (count($parts_msi_ref_line) < 10) trigger_error("MSI ref line with less parts than expected: '".$msi_ref_line."' - any changes in msi ref?", E_USER_ERROR);
			
			//not in filtered bed skip and continue
			if($parts_msi_ref_line[1] != $bed_start || $parts_msi_ref_line[0] != $bed_chr) continue;
			
			//match found: save line and read next lines
			$sites_filtered[] = $msi_ref_line;
			
			//Read next BED line while skipping empty lines
			$bed_filtered_line = "";
			while(! feof($handle_msi_ref_bed_filtered) && nl_trim($bed_filtered_line) == "") 
			{
				$bed_filtered_line = fgets($handle_msi_ref_bed_filtered);
			}
			
			//all lines in filtered file found -> finished
			if (feof($handle_msi_ref_bed_filtered) && nl_trim($bed_filtered_line) == "") break; 
			
			if (count(explode("\t", $bed_filtered_line)) < 3) trigger_error("BED line with less parts than expected: '".$bed_filtered_line."'", E_USER_ERROR);
			list($bed_chr, $bed_start, $bed_end) = explode("\t", $bed_filtered_line);
		}
	}
	
	if (count($sites_filtered) <= 1)  trigger_error("No microsatelite sites overlapping with sample target region!", E_USER_ERROR);
	if (count($sites_filtered) < 500) trigger_error("Very few microsatelite sites overlapping with sample target region! Count:".count($sites_filtered)-1, E_USER_WARNING);	
	
	file_put_contents($tmp_targeted_msi_ref, implode("", $sites_filtered));
	file_put_contents("/mnt/storage2/users/ahott1a1/msi_targeted_test.site", implode("", $sites_filtered));
	$parser->log("finished filtering msi ref file.");
	
	$msi_ref = $tmp_targeted_msi_ref;
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