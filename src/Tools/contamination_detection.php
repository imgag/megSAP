<?php
/** 
	@page contamination_detection
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("contamination_detection", "Detects which sample caused the contamination of a second sample.");
$parser->addString("ps", "Processed sample identifier of contaminated sample.", false);
extract($parser->parse($argv));

//init
$db = DB::getInstance("NGSD");
$ngsbits = get_path("ngs-bits");
$info = get_processed_sample_info($db, $ps);

//default settings for WES/panel
$min_dp = 50;
$r1_start = 0.05;
$r1_end = 0.30;
$r2_start = 0.70;
$r2_end = 0.95;
if ($info['sys_type']=="WGS")
{
	print "##WGS sample > adapting parameters\n";
	$min_dp = 30;
	$r1_start = 0.10;
	$r1_end = 0.20;
	$r2_start = 0.80;
	$r2_end = 0.90;
}
print "##Using parameters: min_dp={$min_dp} range1={$r1_start}-{$r1_end} range2={$r2_start}-{$r2_end}\n";

//get variants in unexpected range
$c_vars = 0;
$vars = array();
$gsvar = substr($info['ps_bam'], 0, -4).".GSvar";
print "##GSvar file: {$gsvar}\n";
$h = fopen($gsvar, "r");
while(!feof($h))
{
	$line = trim(fgets($h));
	if ($line=="" || $line[0]=="#") continue;
	
	list($chr, $start, , $ref, $obs, , , $qual) = explode("\t", $line);
	$qual = explode(";", $qual);
	
	//only autosomes
	if (!is_numeric(substr($chr,3))) continue;
	
	//only SNPs
	if (strlen($ref)!=1 || strlen($obs)!=1 || $ref=="-" || $obs=="-") continue;
	
	//only high-depth variants
	$dp = null;
	foreach($qual as $entry)
	{
		if (starts_with($entry, "DP="))
		{
			$dp = substr($entry, 3);
		}
	}
	if ($dp<$min_dp)
	{
		continue;
	}
	
	++$c_vars;
	
	//only strange AF
	$af = null;
	foreach($qual as $entry)
	{
		if (starts_with($entry, "AF="))
		{
			$af = substr($entry, 3);
		}
	}
	if (($af>$r1_start && $af<$r1_end) || ($af>$r2_start && $af<$r2_end))
	{
		$vars["$chr:$start $ref>$obs"] = true;		
	}
}
fclose($h);
print "##Selected ".count($vars)." variants (".number_format(100.0*count($vars)/$c_vars,2)."% of {$c_vars} autosomal SNPs)\n";

//determine samples of the same processing system
$tmp = temp_file(".tsv");
$sys = $info['sys_name_short'];
$pipeline = [
	["{$ngsbits}NGSDExportSamples", "-system {$sys} -no_bad_runs -no_bad_samples -no_tumor -no_ffpe -add_path"],
	["{$ngsbits}TsvFilter", "-filter 'system_name_short is {$sys}'"], //this is necessary because NGSDExportSamples performs fuzzy match
	["{$ngsbits}TsvSlice", "-cols name,path -out {$tmp}"],
];
$parser->execPipeline($pipeline, "NGSD sample extraction", true);
$samples = file($tmp);
print "##Processing ".count($samples)." samples of the processing system {$sys}...\n";

//process samples
print "#sample\toverlap_perc\thom_perc\n";
foreach($samples as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	list($ps, $path) = explode("\t", $line);
	$gsvar = "{$path}/{$ps}.GSvar";
	
	if (!file_exists($gsvar))
	{
		print "##Skipping sample {$ps}: GSvar file missing {$gsvar}\n";
		continue;
	}
	
	//try to re-find variants in the sample
	$c_refound = 0;
	$c_hom = 0;
	$h = fopen($gsvar, "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || $line[0]=="#") continue;
		
		list($chr, $start, , $ref, $obs, $geno) = explode("\t", $line);
		
		if (isset($vars["$chr:$start $ref>$obs"]))
		{
			++$c_refound;
			
			if ($geno=="hom")
			{
				++$c_hom;
			}
		}
	}
	fclose($h);
	
	print "{$ps}\t".number_format(100.0*$c_refound/count($vars), 2)."\t".number_format(100.0*$c_hom/$c_refound, 2)."\n";
}

?>
