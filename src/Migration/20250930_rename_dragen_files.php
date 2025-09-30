<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("20250930_rename_dragen_files", "Rename DRAGEN output files back to Illumina schema.");
$parser->addString("ps", "Processed sample name.", true, "");
$parser->addString("folder", "Processed sample folder.", true, "");
extract($parser->parse($argv));

//rename DRAGEN folder
if (file_exists("{$folder}/dragen_variant_calls/"))
{
	exec2("mv {$folder}/dragen_variant_calls/ {$folder}/dragen/");
}
else if (file_exists("{$folder}/dragen/"))
{
	print "Nothing to do - DRAGEN folder already renamed!\n";
	exit(0);
}
else
{
	print "Error: DRAGEN folder missing!\n";
	exit(1);
}

//rename files inside DRAGEN folder
list($stdout) = exec2("ls {$folder}/dragen/");
foreach($stdout as $file)
{
	$file = trim($file);
	if ($file=="") continue;
	
	//keep logs folder untouched
	if ($file=="logs") continue;
	
	$filepath = "{$folder}/dragen/{$file}";
	$base = "{$folder}/dragen/{$ps}";
	
	//small variants VCF/gVCF
	if ($file=="{$ps}_dragen.vcf.gz")
	{
		exec2("mv {$filepath} {$base}.hard-filtered.vcf.gz");
	}
	else if ($file=="{$ps}_dragen.vcf.gz.tbi")
	{
		exec2("mv {$filepath} {$base}.hard-filtered.vcf.gz.tbi");
	}
	else if ($file=="{$ps}_dragen.gvcf.gz")
	{
		exec2("mv {$filepath} {$base}.hard-filtered.gvcf.gz");
	}
	else if ($file=="{$ps}_dragen.gvcf.gz.tbi")
	{
		exec2("mv {$filepath} {$base}.hard-filtered.gvcf.gz.tbi");
	}

	//CNVs
	else if ($file=="{$ps}_dragen_cnv.vcf.gz")
	{
		exec2("mv {$filepath} {$base}.cnv.vcf.gz");
	}
	else if ($file=="{$ps}_dragen_cnv.vcf.gz.tbi")
	{
		exec2("mv {$filepath} {$base}.cnv.vcf.gz.tbi");
	}
	else if ($file=="{$ps}_dragen_cnvs.bw")
	{
		exec2("mv {$filepath} {$base}.target.counts.diploid.bw");
	}
	else if ($file=="{$ps}_dragen_cnvs.vcf.gz")
	{
		exec2("mv {$filepath} {$base}.cnv.vcf.gz");
	}
	else if ($file=="{$ps}_dragen_cnvs.vcf.gz.tbi")
	{
		exec2("mv {$filepath} {$base}.cnv.vcf.gz.tbi");
	}
	
	//SVs
	else if ($file=="{$ps}_dragen_svs.vcf.gz")
	{
		exec2("mv {$filepath} {$base}.sv.vcf.gz");
	}
	else if ($file=="{$ps}_dragen_svs.vcf.gz.tbi")
	{
		exec2("mv {$filepath} {$base}.sv.vcf.gz.tbi");
	}
	
	//unexpected file > abort
	else
	{
		print "UNEXPECTED: {$filepath}\n";
		die(1);
	}
}

?>