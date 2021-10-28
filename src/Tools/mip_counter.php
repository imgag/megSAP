<?php

/** 
	@page mip_counter
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("mip_counter", "Extracts MIP counts from a sample (on-target and off-target).");
$parser->addString("ps", "Processed sample name.", false);
//optional
$parser->addInfile("add_bed", "Additional regions BED file to add to target region, e.g. for QC regions.", true);
extract($parser->parse($argv));

function get_arm_counts($bam, $roi, $fastq_r1, $arm_len)
{
	global $parser;
	
	$output = array();
	
	//extract reads IDs
	$read_ids = array();
	list($tmp) = $parser->exec(get_path("samtools")." view", "-L $roi -M $bam", true);
	foreach($tmp as $line)
	{
		$parts = explode(":", $line);
		$id = implode(":", array_slice($parts, 0, 7)); //first seven parts are relevant for Illumina reads
		$read_ids[] = $id;
	}
	$ids_file = $parser->tempFile(".txt");
	file_put_contents($ids_file, join("\n", $read_ids));
	
	//extract reads with IDs
	$fastq_out = $parser->tempFile(".txt");
	list($read_data) = $parser->exec(get_path("ngs-bits")."FastqExtract", "-in {$fastq_r1} -ids {$ids_file} -out {$fastq_out}", true);
	$i=0;
	$h = gzopen2($fastq_out, "r");
	while(!feof($h))
	{
		$line = fgets($h);
		if ($i==1)
		{
			$arm = substr($line, 0, $arm_len);
			@$output[$arm] += 1;
		}

		++$i;
		if ($i==4) $i=0;
	}
	
	arsort($output);
	
	return $output;
}

function arm_to_name($arm, $arm2name)
{
	foreach($arm2name as $curr_arm => $curr_key)
	{
		if (contains($curr_arm, $arm)) return $curr_key;
	}
	
	$arm = rev_comp($arm);
	foreach($arm2name as $curr_arm => $curr_key)
	{
		if (contains($curr_arm, $arm)) return $curr_key." (reverse-complementary)";
	}
	
	return "n/a";
}

//init
$db = DB::getInstance("NGSD");
$data_folder = get_path("data_folder");
$info = get_processed_sample_info($db, $ps);
$bam = $info['ps_bam'];

//construct target region
$sys = load_system($system, $ps);
$roi = $sys['target_file'];
if (isset($add_bed))
{
	$tmp = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedAdd", "-in {$roi} {$add_bed} -out {$tmp}", true);
	$roi = $tmp;
}
$tmp = $parser->tempFile(".bed");
$parser->exec(get_path("ngs-bits")."BedExtend", "-in {$roi} -n 200 -out {$tmp}", true);
$roi = $tmp;

//contruct off-target region
$not_roi = $parser->tempFile(".bed");
$parser->exec(get_path("ngs-bits")."BedSubtract", "-in {$data_folder}/enrichment/WGS_hg38.bed -in2 {$roi} -out {$not_roi}", true);

//determin shortest MIP extension arm length
$arm_len = 999;
$arm2name = array();
$mip_file = get_path("data_folder")."/mipfiles/".$sys["name_short"].".txt";
$tmp = file($mip_file);
foreach($tmp as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]==">") continue;
	$parts = explode("\t", $line);
	$arm = $parts[10];
	$arm_len = min($arm_len, strlen($arm));
	
	$name = $parts[19];
	if (isset($arm2name[$arm])) trigger_error("MIP arm '$arm' is not unique - it matches MIP '$name' and '".$arm2name[$arm]."'!", E_USER_WARNING);
	$arm2name[$arm] = $name;
}

//determine FASTQ file of read 1
$ps_folder = $info['ps_folder'];
$tmp = glob("{$ps_folder}/*R1*.fastq.gz");
$fastq_r1 = $tmp[0];

//print on-target arms
print "on-target MIPs:\n";
$arms = get_arm_counts($bam, $roi, $fastq_r1, $arm_len);
foreach($arms as $arm => $count)
{
	print "{$count}\t{$arm}\t".arm_to_name($arm, $arm2name)."\n";
}

//print off-target arms
print "\n";
print "off-target MIPs:\n";
$arms = get_arm_counts($bam, $not_roi, $fastq_r1, $arm_len);
foreach($arms as $arm => $count)
{
	print "{$count}\t{$arm}\t".arm_to_name($arm, $arm2name)."\n";
}
?>