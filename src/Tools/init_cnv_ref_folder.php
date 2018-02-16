<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

/**
	@page init_cnv_ref_folder
*/

// parse command line arguments
$parser = new ToolBase("init_cnv_ref_folder", "Initialize coverage folder of a processing system with data from all valid samples.");
$parser->addString("name", "Processing system short name.", false);
$parser->addFlag("clear", "Remove existing coverage files.");
extract($parser->parse($argv));

//check system
$db = DB::getInstance("NGSD");
$res = $db->executeQuery("SELECT * FROM processing_system WHERE name_short=:sys", array("sys"=>$name));
if (count($res)!=1)
{
	trigger_error("Invalid processing system short name '$name'! Valid names are: ".implode(", ", $db->getEnum("processing_system", "name_short")), E_USER_ERROR);	
}
$roi = trim($res[0]['target_file']);
if ($roi=="")
{
	trigger_error("Processing system '$name' has no target region file annotated!", E_USER_ERROR);	
}

//create/clear folder
$ref_folder = get_path("data_folder")."/coverage/$name/";		
if (!is_dir($ref_folder))
{
	mkdir($ref_folder);
	if (!chmod($ref_folder, 0777))
	{
		trigger_error("Could not change privileges of folder '{$ref_folder}'!", E_USER_ERROR);
	}
}
if ($clear)
{
	exec2("rm -rf $ref_folder/*.cov");
}

//print output
print "Processing system: $name\n";
print "Target file: $name\n";
print "Coverage folder: $ref_folder\n";
print "\n";

//check samples
list($samples) = exec2(get_path("ngs-bits")."NGSDExportSamples -sys {$name} -quality bad -check_path | cut -f1,18 | grep -v ps.name");
$bams = array();
$valid = 0;
foreach($samples as $line)
{
	list($sample, $path) = explode("\t", trim($line));
	//check sample is valid
	if(!is_valid_ref_sample_for_cnv_analysis($sample)) continue;
	++$valid;
	//check bam
	if (!starts_with($path, "found - ")) continue;
	$bam = substr($path, 8)."/".$sample.".bam";
	if (file_exists($bam))
	{
		$bams[] = $bam;
	}
}
print "Found ".count($samples)." samples for this processing system in NGSD.\n";
print "Found $valid samples that are valid reference samples.\n";
print "Found ".count($bams)." samples with BAM files.\n";

//create new coverage files
foreach($bams as $bam)
{
	$cov_file = "$ref_folder/".basename($bam, ".bam").".cov";
	
	//skip existing coverage files
	if (file_exists($cov_file))
	{
		print "$bam skipped (coverage file already exists)\n";
		continue;
	}
	
	//calculate coverage file
	print "$bam processing...\n";
	exec2(get_path("ngs-bits")."BedCoverage -bam $bam -in $roi -out $cov_file -min_mapq 0");
}

//chmod
exec2("chmod 775 $ref_folder/*.cov");

?>
