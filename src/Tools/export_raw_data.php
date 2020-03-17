<?php
/** 
	@page export_raw_data 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_raw_data", "Exports RAW data.");
$parser->addStringArray("samples", "Processed sample names.", false);
$parser->addString("out", "Output folder name.", false);
extract($parser->parse($argv));

//init
$db = DB::getInstance("NGSD");
$ngsbits = get_path("ngs-bits");

//check folder
if (!file_exists($out) || !is_dir($out))
{
	print "Output folder '$out' is missing - creating it ...\n";
	mkdir($out, 0777, true);
}

//check and create meta data file
print "Checking samples ...\n";
$meta = [];
$meta[] = "#sample\texternal_name\tgender\tdisease_group_and_status\tkit\tsequencer\n";
foreach($samples as $ps)
{
	print "  $ps\n";

	$info = get_processed_sample_info($db, $ps);
	$quality = $info['ps_quality'];
	if ($quality=="bad") trigger_error("Sample '$ps' has 'bad' quality!", E_USER_ERROR);
	$bam = $info['ps_bam'];
	if (!file_exists($bam)) trigger_error("Sample '$ps' BAM file is missing: $bam", E_USER_ERROR);
	
	$gender = $info['gender'];
	$external = $info['name_external'];
	$device = $info['device_type'];
	$system = $info['sys_name'];
	$d_group = $info['disease_group'];
	$d_status = $info['disease_status'];
	$meta[] = "{$ps}\t{$external}\t{$gender}\t{$d_group} ({$d_status})\t{$system}\t{$device}\n";
}
file_put_contents("$out/meta_data.tsv", $meta);

//export FASTQ data
print "Creating FASTQ files ...\n";
foreach($samples as $ps)
{
	print "  $ps\n";
	$info = get_processed_sample_info($db, $ps);
	$bam = $info['ps_bam'];
	
	exec2("{$ngsbits}/BamToFastq -in {$bam} -out1 {$out}/{$ps}_R1_001.fastq.gz -out2 {$out}/{$ps}_R2_001.fastq.gz");
}

//zip
print "Zipping output folder...\n";
exec2("zip -r {$out}.zip $out"); //TODO -e for excryption?

//move
print "You can move the file to the webserver using:\n";
$user = exec('whoami');
print "scp {$out}.zip {$user}@imgag.de:/var/www/html/download/{$user}/\n";


?>