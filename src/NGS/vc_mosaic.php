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
$parser->addEnum("type", "Processing system type: WGS or WES", false, ["WES", "WGS"]);

//optional
$parser->addInt("extend", "Extend gene regions by 'extend' bases.", true, 20);
//filter parameter:
$parser->addFloat("max_af", "Maximum allele-frequency of a variant to be considered as a mosaic", true, 0.5);
$parser->addInt("min_obs", "Minimum observation per strand. if not given is decided by the type parameter: WGS = 1; WES = 2", true, -1);
$parser->addFloat("max_gnomad_af", "Maximum allowed allel frequency in gnomad", true, 0.01);
//freebayes calling parameters:
$parser->addFloat("min_af", "Minimum allele-frequency", true, 0.01);
$parser->addInt("min_mq", "Minimum mapping quality (freebayes)", true, 55);
$parser->addInt("min_bq", "Minimum base quality (freebayes)", true, 25);
$parser->addInt("min_qsum", "Minimum alternate quality sum (freebayes)", true, 90);

$parser->addInt("threads", "Number of threads in annotation step", true, 1);

//debug parameters
$parser->addFlag("debug", "Enable additional output for debugging.");
extract($parser->parse($argv));

// Validate input:
if ($max_af <= 0)
{
	trigger_error("The maximum allel frequency cannot be smaller or equal to zero.", E_USER_ERROR);
}

if ($gnomad_af <= 0)
{
	trigger_error("The maximum allowed gnomad allel frequency cannot be smaller or equal to zero.", E_USER_ERROR);
}

if ($min_obs == -1)
{
	if ($type == "WES")
	{
		$min_obs=2;
	}
	else if ($type == "WGS")
	{
		$min_obs=1;
	}
	else
	{
		trigger_error("Unknown type: $type , cannot automatically determine appropriate min_obs parameter.", E_USER_ERROR);
	}
}


/*
Filters the base results and removes called variants that:
	- already occur in the given .vcf file from the sample.
	- have a higher allel frequency than max_af.
	- occure on one strand less than min_obs.
	- have a higher allel frequency in gnomAD than max_gnomad_af.
	- occure in low confidence regions.
*/

function filter_vcf($vcf, $out, $sample_vcf, $max_af, $min_obs, $min_quality, $max_gnomad_af)
{
	$result = [];
	
	$vars_germline = array();

	$h = gzopen($sample_vcf, "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || $line[0]=="#") continue;
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $all_info, $format, $format_values) = explode("\t", $line);	
		$mut = "$chr:$pos $ref>$alt";
		$vars_germline[$mut] = true;
	}
	fclose($h);
	
	
	foreach (file($vcf) as $line)
	{
		$line = trim($line);
		
		if (starts_with($line, "#"))
		{
			$result[] = $line;
			continue;
		}
		

		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $all_info, $format, $format_values) = explode("\t", $line);	
		$var = "$chr:$pos $ref>$alt";

		list($saf, $sar, $af, $gnomad_af) = parse_vcf_info($all_info, $format, $format_values);
		
		if ($af > $max_af)
		{
			continue;
		}
		
		if ($saf < $min_obs || $sar < $min_obs)
		{
			continue;
		}
		
		if (floatval($gnomad_af) > $max_gnomad_af)
		{
			continue;
		}
		
		if (str_contains($filter, "low_conf_region"))
		{
			continue;
		}
		
		if (array_key_exists($var, $vars_germline))
		{
			continue;
		}
		
		$result[] = $line;
	}
	
	file_put_contents($out, implode("\n", $result));
}

function parse_vcf_info($info, $format, $format_values)
{
	$saf = "";	// alt occurences on forward strand
	$sar = "";	// alt occurences on reverse strand
	$gnomad_af = "";
	
	$info = explode(";", $info);
	foreach($info as $i)
	{
		if(starts_with($i, "SAF="))
		{
			$saf = substr($i, 4);
		}
		if(starts_with($i, "SAR="))
		{
			$sar = substr($i, 4);
		}
		if(starts_with($i, "gnomADg_AF="))
		{
			$gnomad_af = substr($i, 11);
		}
	}
	
	$af = (floatval($sar) + floatval($saf)) / floatval($dp);
	return [$saf, $sar, $af, $gnomad_af];
}


//**MAIN**//

$time_start = microtime(true);

//create target region
$gene_regions = temp_file("_gene_regions.bed");
$parser->exec(get_path("ngs-bits")."/GenesToBed", "-in $genes -source ccds -mode exon -fallback -out $gene_regions", "Creating target region");
if ($debug) print "Target regions: $gene_regions \n";

$final_region = $gene_regions;
//extend target region
if ($extend > 0)
{
	$bed_ext = temp_file("_extended.bed");
	$parser->exec(get_path("ngs-bits")."/BedExtend", "-in $gene_regions -n $extend -out $bed_ext", "Extending target region"); 
	if ($debug) print "Regions extended: $bed_ext \n";


	$bed_merged = temp_file("_merged.bed");
	$parser->exec(get_path("ngs-bits")."/BedMerge", "-in $bed_ext -out $bed_merged", "Merging target region"); 
	if ($debug) print "Regions merged: $bed_merged \n";
	$final_region = $bed_merged;
}

//call variants
$called_vcf = temp_file(".vcf");
exec(get_path("freebayes")." -t $final_region -b $in -f /tmp/local_ngs_data//GRCh38.fa --pooled-continuous --min-alternate-fraction $min_af --min-mapping-quality $min_mq --min-base-quality $min_bq --min-alternate-qsum $min_qsum > $out");
if ($debug) print "variant calls: $called_vcf\n";


//**POST PROCESSING**//

//split complex variants to primitives
$tmp_split_vars = temp_file("_split_vars.vcf");
exec("cat $vcf | ".get_path("vcflib")."vcfallelicprimitives -kg > $tmp_split_vars");
if ($debug) print "Split complex variants: $tmp_split_vars\n";

// split multi-allelic variants
$tmp_break_multi = temp_file("_break_multi.vcf");
exec("cat $tmp_split_vars | ".get_path("vcflib")."vcfbreakmulti > $tmp_break_multi");
if ($debug) print "Split multi-allelic variants: $tmp_break_multi\n";

// normalize all variants and align INDELs to the left
$tmp_left_norm = temp_file("_left_norm.vcf");
exec(get_path("ngs_bits")."VcfLeftNormalize -in $tmp_break_multi -out $tmp_left_norm");
if ($debug) print "Left normalized variants: $tmp_left_norm\n";

// //annotate the vcf file
$tmp_annotated = temp_file("_annotated.vcf");
exec(get_path("ngs_bits")."VcfAnnotateFromVcf -in $tmp_left_norm -out $tmp_annotated -annotation_file /mnt/storage2/GRCh38/share/data/dbs/gnomAD/gnomAD_genome_v3.1.1_GRCh38.vcf.gz -info_ids AF -id_prefix gnomADg -threads $threads");
if ($debug) print "Annotated variants: $tmp_annotated\n";

// filter: sar saf (min observations)
$tmp_filtered = temp_file("_filtered.vcf");
filter_vcf($tmp_annotated, $tmp_filtered, $vcf, $max_af, $min_obs, $min_quality, $max_gnomad_af);
if ($debug) print "Filtered variants: $tmp_filtered\n";

//sort variants by genomic position
$tmp_sorted_vcf = temp_file("_sorted.vcf");
exec(get_path("ngs_bits")."VcfSort -in $tmp_filtered -out $tmp_sorted_vcf");
if ($debug) print "Sorted variants: $tmp_sorted_vcf\n";

// fix error in VCF file and strip unneeded information
exec("cat $tmp_sorted_vcf | php ".execTool("NGS/vcf_fix.php", "--keep_wt_calls > $out");

$time_end = microtime(true);
$time = $time_end - $time_start;
if ($debug) print "Done - total time taken: $time s \n";

?>