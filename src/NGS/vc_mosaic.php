<?php

/**
  @page vc_mosaic
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_mosaic", "Call mosaic mutations in bam file. Creates a VCF.gz file.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addInfile("vcf", "VCF(.gz) with mutations called following the diploid model", false);
$parser->addOutfile("out", "Output VCF file.", false);
$parser->addEnum("type", "Processing system type: WGS or WES", false, ["WES", "WGS"]);

//optional
$parser->addInfile("genes", "File with gene names in whose exonic regions the mutations are searched for.", true, "");
$parser->addInt("extend", "Extend gene regions by 'extend' bases.", true, 20);
$parser->addInfile("target", "File with genome regions in which the mutations are searched for. (will not use extend parameter)", true, "");
$parser->addString("build", "The genome build to use.", true, "GRCh38");
//filter parameter:
$parser->addFloat("max_af", "Maximum allele-frequency of a variant to be considered as a mosaic", true, 0.5);
$parser->addInt("min_obs", "Minimum observation per strand. if not given is decided by the type parameter: WGS = 1; WES = 2", true, -1);
$parser->addFloat("max_gnomad_af", "Maximum allowed allel frequency in gnomad", true, 0.01);
//freebayes calling parameters:
$parser->addFloat("min_af", "Minimum allele-frequency (freebayes)", true, 0.01);
$parser->addInt("min_mq", "Minimum mapping quality (freebayes)", true, 50);
$parser->addInt("min_bq", "Minimum base quality (freebayes)", true, 25);
$parser->addInt("min_qsum", "Minimum alternate quality sum (freebayes)", true, 90);

$parser->addInt("threads", "Number of threads used.", true, 1);

//debug parameters
$parser->addFlag("debug", "Enable additional output for debugging.");
$parser->addInfile("called_vcf", "A vcf that was already called for mosaics and only post-processing should be done.", true, "");
extract($parser->parse($argv));

$time_start = microtime(true);

if (($genes == "" && $target == "") || ($genes != "" && $target != ""))
{
	trigger_error("vc_mosaic: Needs either gene parameter or target parameter to be set. Don't give both or neither.", E_USER_ERROR);
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

// Validate input:
if ($max_af <= 0)
{
	trigger_error("The maximum allel frequency cannot be smaller or equal to zero.", E_USER_ERROR);
}

if ($max_gnomad_af <= 0)
{
	trigger_error("The maximum allowed gnomad allel frequency cannot be smaller or equal to zero.", E_USER_ERROR);
}

if ($min_obs < 0)
{
	trigger_error("The minimum observations cannot be negative.", E_USER_ERROR);
}

/*
Filters the base results and removes called variants that:
	- already occur in the given .vcf file from the sample.
	- have a higher allel frequency than max_af.
	- occure on one strand less than min_obs.
	- have a higher allel frequency in gnomAD than max_gnomad_af.
	- occure in low confidence regions.
*/

function filter_vcf($vcf, $out, $sample_vcf, $max_af, $min_obs, $max_gnomad_af)
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
		
		if (contains($filter, "low_conf_region"))
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
	$gnomad_af = "0";
	$dp = "";
	
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
		if(starts_with($i, "DP="))
		{
			$dp = substr($i, 3);
		}
		if(starts_with($i, "gnomADg_AF="))
		{
			$gnomad_af = substr($i, 11);
		}
	}
	
	if ($saf == "" || $sar == "")
	{
		
		var_dump($info);
		trigger_error("Couldn't parse needed info values (at least one of SAR, SAF)", E_USER_ERROR);
	}
	
	if ($dp == "")
	{
		//parse depth from the format values:
		$dp_idx = -1;
		$format_ids = explode(":", $format);
		for($i=0; $i<count($format_ids); $i++)
		{
			$id = $format_ids[$i];
			// echo "$i . ID: $id\n";
			if ($id == "DP")
			{
				// echo ""
				$dp_idx = $i;
			}
		}
		if ($dp_idx == -1)
		{
			var_dump($format_ids);
			trigger_error("No depth value found in info fields and not annotated in the format values in the vcf: format - $format .", E_USER_ERROR);
		}
		$dp = explode(":", $format_values)[$dp_idx];
	}
	
	$saf = intval($saf);
	$sar = intval($sar);
	$dp = intval($dp);
	$gnomad_af = floatval($gnomad_af);
	
	
	$af = ($sar + $saf) / $dp;
	

	return [$saf, $sar, $af, $gnomad_af];
}


if ($called_vcf == "")
{

//**MAIN**//
if ($debug) print "starting main\n";

$final_region = "";
if ($genes != "")
{
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
}
else
{
	$final_region = $target;
}
$freebayes_start = microtime(true);
//call variants
$called_vcf = temp_file(".vcf");
$parser->exec("php ".repository_basedir()."src/NGS/vc_freebayes.php ", " -target $final_region -bam $in -out $called_vcf -build $build -no_ploidy -min_af $min_af -min_mq $min_mq -min_bq $min_bq -min_qsum $min_qsum -raw_output");

$freebayes_end = microtime(true);
if ($debug) print "variant calling took ".($freebayes_end-$freebayes_start)."s: $called_vcf\n";
}

//**POST PROCESSING**//

$post_processing_start = microtime(true);
//split complex variants to primitives
$tmp_split_vars = temp_file("_split_vars.vcf");
exec("cat $called_vcf | ".get_path("vcflib")."vcfallelicprimitives -kg > $tmp_split_vars");
if ($debug) print "Split complex variants: $tmp_split_vars\n";

// split multi-allelic variants
$tmp_break_multi = temp_file("_break_multi.vcf");
exec("cat $tmp_split_vars | ".get_path("vcflib")."vcfbreakmulti > $tmp_break_multi");
if ($debug) print "Split multi-allelic variants: $tmp_break_multi\n";

// normalize all variants and align INDELs to the left
$tmp_left_norm = temp_file("_left_norm.vcf");
exec(get_path("ngs-bits")."VcfLeftNormalize -in $tmp_break_multi -out $tmp_left_norm");
if ($debug) print "Left normalized variants: $tmp_left_norm\n";

// //annotate the vcf file
$tmp_annotated = temp_file("_annotated.vcf");
exec(get_path("ngs-bits")."VcfAnnotateFromVcf -in $tmp_left_norm -out $tmp_annotated -annotation_file /mnt/storage2/GRCh38/share/data/dbs/gnomAD/gnomAD_genome_v3.1.1_GRCh38.vcf.gz -info_ids AF -id_prefix gnomADg -threads $threads");
if ($debug) print "Annotated variants: $tmp_annotated\n";

// filter
$tmp_filtered = temp_file("_filtered.vcf");
filter_vcf($tmp_annotated, $tmp_filtered, $vcf, $max_af, $min_obs, $max_gnomad_af);
if ($debug) print "Filtered variants: $tmp_filtered\n";

//sort variants by genomic position
$tmp_sorted_vcf = temp_file("_sorted.vcf");
exec(get_path("ngs-bits")."VcfSort -in $tmp_filtered -out $tmp_sorted_vcf");
if ($debug) print "Sorted variants: $tmp_sorted_vcf\n";

// fix error in VCF file and strip unneeded information
$tmp_fixed_vcf = temp_file("_fixed.vcf");
exec("cat $tmp_sorted_vcf | php ".repository_basedir()."src/NGS/vcf_fix.php --mosaic_mode | bgzip > $out");

$time_end = microtime(true);
if ($debug) print "Post processing done - time taken: ".($time_end-$post_processing_start)." s \n";
if ($debug) print "Done - total time taken: ".($time_end - $time_start)." s \n";

?>