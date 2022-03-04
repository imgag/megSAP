<?php

/**
  @page mosaic_finder
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mosaic_finder", "Call mosaic mutations in bam file. Creates a VCF file.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addInfile("vcf", "VCF(.gz) with mutations called following the diploid model", false);
$parser->addInfile("genes", "File with gene names in whose regions the mutations are searched for.", false);
$parser->addOutfile("out", "Output VCF.gz file.", false);

//optional
$parser->addFloat("min_af", "Minimum allele-frequency", true, 0.01);
$parser->addFloat("max_af", "Maximum allele-frequency", true, 1.0);
$parser->addInt("min_obs", "Minimum observation per strand.", true, 2);

//debug parameters
$parser->addFlag("debug", "Enable additional output for debugging.");
$parser->addString("debug_folder", "Save temporary files in given folder instead of /tmp/.", true, "");
extract($parser->parse($argv));

if ($debug)
{
	echo "Input variables:\n";
	echo $in."\n";
	echo $vcf."\n";
	echo $genes."\n";
	echo $out."\n";
	echo $min_af."\n";
	echo $max_af."\n";
	echo $min_obs."\n";
}

if ($debug_folder != "")
{
	if ( ! is_dir($debug_folder))
	{
		trigger_error("Provided debug folder doesn't exist or is not a directory.", E_USER_ERROR);
	}
}

//TODO: input testing - everything ok to start?

//create target region
$gene_regions = temp_file("_gene_regions.bed");
if ($debug_folder != "")
{
	$gene_regions = $debug_folder."/gene_regions.bed";
}

$parser->exec(get_path("ngs-bits")."/GenesToBed", "-in $genes -source ccds -mode exon -fallback -out $gene_regions", "Creating target region");
if ($debug) print "Target regions: $gene_regions \n";


//extend target region
$bed_ext = temp_file("_extended.bed");
if ($debug_folder != "")
{
	$bed_ext = $debug_folder."/extended_regions.bed";
}
$parser->exec(get_path("ngs-bits")."/BedExtend", "-in $gene_regions -n 20 -out $bed_ext", "Extending target region"); 
if ($debug) print "Regions extended: $bed_ext \n";


$bed_merged = temp_file("_merged.bed");
if ($debug_folder != "")
{
	$bed_merged = $debug_folder."/merged_regions.bed";
}
$parser->exec(get_path("ngs-bits")."/BedMerge", "-in $bed_ext -out $bed_merged", "Merging target region"); 
if ($debug) print "Regions merged: $bed_merged \n";

//load germline variants
$vars_germline = array();

$h = gzopen($vcf, "r");
while(!feof($h))
{
	$line = trim(fgets($h));
	if ($line=="" || $line[0]=="#") continue;
	list($chr, $pos, , $ref, $alt) = explode("\t", $line, 6);
	$vars_germline["$chr:$pos $ref>$alt"] = true;
}
fclose($h);
if ($debug) print "Vcf read. Count ".count($vars_germline)."\n";


//call variants
$called_vcf = temp_file(".vcf");
if ($debug_folder != "")
{
	$called_vcf = $debug_folder."/called_vars.vcf";
}
exec(get_path("freebayes")." -t $bed_merged -b $in -f /tmp/local_ngs_data//GRCh38.fa --pooled-continuous --min-alternate-fraction $min_af --min-mapping-quality 50 --min-base-quality 25 --min-alternate-qsum 90 > $called_vcf");
if ($debug) print "variant calls: $called_vcf\n";

//filter variants
$output = array();
$file = file($called_vcf);
foreach($file as $line)
{
	if (starts_with($line, "#"))
	{
		$output[] = trim($line);
		continue;
	}
	
	list($chr, $pos, , $ref, $alt, $qual, , $info, $format, $format_values) = explode("\t", $line);
	
	$saf = "";
	$sar = "";
	$dp = "";
	$info = explode(";", $info);
	foreach($info as $i)
	{
		if(starts_with($i, "DP="))
		{
			$dp = substr($i, 3);
		}
		if(starts_with($i, "SAF="))
		{
			$saf = substr($i, 4);
		}
		if(starts_with($i, "SAR="))
		{
			$sar = substr($i, 4);
		}
	}
	$alt = explode(",", $alt);
	$saf = explode(",", $saf);
	$sar = explode(",", $sar);
	
	for($i=0; $i<count($alt); ++$i)
	{
		// skip variants that were called by normal germline variant calling
		$tag = "$chr:$pos $ref>".$alt[$i];
		if (isset($vars_germline[$tag]))
		{
			continue;
		}
		//at least n observeration in each direction
		if ($saf[$i]<$min_obs || $sar[$i]<$min_obs)
		{
			continue;
		}
		$af = (floatval($saf[$i]) + floatval($sar[$i])) / floatval($dp); 
		if ($af < $min_af || $af > $max_af)
		{
			continue;
		}
		
		$output[] = $line;
		break;
	}
}
$vcf_filtered = temp_file("_filtered.vcf");
if ($debug_folder != "")
{
	$vcf_filtered = $debug_folder."/filtered_vars.vcf";
}
file_put_contents($vcf_filtered, implode("\n", $output));


//normalize variants
$vcf_norm = temp_file("_normalized.vcf");
if ($debug_folder != "")
{
	$vcf_norm = $debug_folder."/normalized_vars.vcf";
}
$parser->exec(get_path("ngs-bits")."/VcfLeftNormalize", "-in $vcf_filtered -out $vcf_norm", "Normalizing variants"); 
if ($debug) print "variant calls - normalized: $vcf_norm \n";

//annotate variants
$parser->execTool("NGS/an_vep.php", "-in $vcf_norm -out $out -ps_name ".basename($in, ".bam")." -no_splice", "Annotating variants with VEP", "NGS");






?>