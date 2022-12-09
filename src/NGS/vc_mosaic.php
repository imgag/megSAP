<?php

/**
  @page vc_mosaic
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_mosaic", "Call mosaic variants in a BAM file. Creates a VCF.GZ file.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addOutfile("out", "Output VCF file.", false);
$parser->addInt("min_obs", "Minimum observation per strand. Recommended value for WGS = 1 and WES = 2", false);
$parser->addInfile("target", "File with target region for variant calling.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFloat("max_af", "Maximum allele-frequency of a variant to be considered as a mosaic", true, 0.5);
$parser->addFloat("max_gnomad_af", "Maximum allowed population allele frequency in gnomAD", true, 0.01);
$parser->addFloat("min_af", "Minimum allele frequency (freebayes)", true, 0.01);
$parser->addInt("min_mq", "Minimum mapping quality (freebayes)", true, 50);
$parser->addInt("min_bq", "Minimum base quality (freebayes)", true, 25);
$parser->addInt("min_qsum", "Minimum alternate quality sum (freebayes)", true, 90);
$parser->addInt("threads", "Number of threads used.", true, 1);
$parser->addFlag("no_zip", "Do not gzip output VCF file.");
extract($parser->parse($argv));

//check parameters
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

//init
$genome = genome_fasta($build);

function filter_vcf($vcf, $out, $max_af, $min_obs, $max_gnomad_af)
{
	$result = [];
	
	//filter mosaic calls:
	foreach (file($vcf) as $line)
	{
		$line = trim($line);
		
		if (starts_with($line, "#"))
		{
			$result[] = $line;
			continue;
		}
		
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $all_info, $format, $format_values) = explode("\t", $line);	
		list($saf, $sar, $af, $gnomad_af) = parse_vcf_info($all_info, $format, $format_values);
		
		if ($af > $max_af)
		{
			continue;
		}

		if ($saf < $min_obs || $sar < $min_obs)
		{
			continue;
		}
		
		if ($gnomad_af > $max_gnomad_af)
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


//call variants
$freebayes_start = microtime(true);
$called_vcf = temp_file(".vcf");
$parser->exec("php ".repository_basedir()."src/NGS/vc_freebayes.php ", " -target $target -bam $in -out $called_vcf -build $build -no_ploidy -min_af $min_af -min_mq $min_mq -min_bq $min_bq -min_qsum $min_qsum -raw_output -threads $threads");

//normalization and annotation
$pipeline = [];
$pipeline[] = array("cat", "$called_vcf");
$pipeline[] = array(get_path("vcflib")."vcfallelicprimitives", "-kg");
$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");
$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize", "-stream -ref $genome");
$tmp_annotated = temp_file("_annotated.vcf");
$gnomad_file = get_path("data_folder")."/dbs/gnomAD/gnomAD_genome_v3.1.2_GRCh38.vcf.gz";
$pipeline[] = array(get_path("ngs-bits")."VcfAnnotateFromVcf", "-out $tmp_annotated -source $gnomad_file -info_keys AF -prefix gnomADg -threads $threads");
$parser->execPipeline($pipeline, "vc_mosaic post processing");

//filter
$tmp_filtered = temp_file("_filtered.vcf");
filter_vcf($tmp_annotated, $tmp_filtered, $max_af, $min_obs, $max_gnomad_af);

//sort variants by genomic position
$tmp_sorted_vcf = temp_file("_sorted.vcf");
$parser->exec(get_path("ngs-bits")."VcfSort", "-in $tmp_filtered -out $tmp_sorted_vcf");

// fix error in VCF file and strip unneeded information
$tmp_fixed_vcf = temp_file("_fixed.vcf");
exec("cat $tmp_sorted_vcf | php ".repository_basedir()."src/NGS/vcf_fix.php --mosaic_mode ".($no_zip ? "" : "| bgzip")." > $out");

?>
