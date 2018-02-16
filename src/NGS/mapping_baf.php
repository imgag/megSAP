<?php

/** 
	@page mapping_baf
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("mapping_baf", "Generate SEG file that contains bafs.");
$parser->addInfile("in",  "Input BAM file.", false);
$parser->addOutfile("out",  "Output IGV file.", false);
//optional
$parser->addInfile("target",  "Target bed file used to identify common SNPs in target region. If unspecified settings for WGS will be used. If -n_in is specified same target file will be used for this sample.", true);
$parser->addInfile("n_in",  "Input normal file in bam format (somatic).", true);
$parser->addFloat("min_af", "Minimum allele frequency of SNPs to use (in 1000g data).", true, 0.01);
$parser->addInt("min_dp", "Minimum depth of SNPs in BAM.", true, 20);
$parser->addFlag("full_wgs", "Do not downsample variant list for WGS.", true);
$parser->addString("base_vcf", "Base variant list, must contain AF field.", true, get_path("data_folder")."dbs/1000G/1000g_v5b.vcf.gz");
extract($parser->parse($argv));

// check if base VCF exists
if (!file_exists($base_vcf))
{
	trigger_error("Base VCF '{$base_vcf}' does not exist!", E_USER_ERROR);
}

$is_somatic = isset($n_in);
$ps_name = basename($in, ".bam");
if ($is_somatic)
{
	$nor_name = basename($n_in, ".bam");
}

// filtered variants
// variant list depends on target region and specified minimum allele fraction
if (isset($target))
{
	$target_shortname = basename($target, ".bed");
}
else if ($full_wgs)
{
	$target_shortname = "WGS_full";
}
else
{
	$target_shortname = "WGS";
}
$filtered_variants = sprintf("%s/baf_%s_%s_af%s.tsv",
	dirname($base_vcf),
	basename($base_vcf, ".vcf.gz"),
	$target_shortname,
	$min_af);

if (!file_exists($filtered_variants))
{
	trigger_error("Filtered variant list does not exist, creating at '{$filtered_variants}'.");

	// filter base VCF by target region (only for enrichment)
	if (isset($target))
	{
		$tmp_filtered_by_region = $parser->tempFile("_filtered_by_region.vcf");
		$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in {$base_vcf} -reg {$target} -out {$tmp_filtered_by_region}", true);
	}
	else
	{
		$tmp_filtered_by_region = $base_vcf;
	}

	// filter known variants (SNPs with high population AF)
	$count = 0;

	$handle = gzopen($tmp_filtered_by_region, "r");
	if ($handle === FALSE)
	{
		trigger_error("Could not open file $known_variants for reading.", E_USER_ERROR);
	}

	$handle_out = fopen($filtered_variants, "w");
	if ($handle_out === FALSE)
	{
		trigger_error("Could not open file $filtered_variants for writing.", E_USER_ERROR);
	}
	fwrite($handle_out, "#chr\tstart\tend\tref\tobs\tid\n");

	while (!feof($handle))
	{
		$line = nl_trim(fgets($handle));

		// skip empty lines and comments
		if ($line == "" || starts_with($line, "#"))
		{
			continue;
		}

		$cols = explode("\t", $line);
		if (count($cols) < 8)
		{
			trigger_error("VCF file line contains less than 8 columns:\n$line", E_USER_ERROR);
		}

		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = $cols;

		// skip special chromosomes
		if (!chr_check($chr, 22, false))
		{
			continue;
		}

		// skip indels
		if (strlen($ref) > 1 || strlen($alt) > 1)
		{
			continue;
		}

		// skip low AF variants
		$af = 0.0;
		foreach (explode(";", $info) as $entry)
		{
			if (starts_with($entry, "AF="))
			{
				$af = substr($entry, 3);
				break;
			}
		}
		if ($af < $min_af)
		{
			continue;
		}

		// WGS: use only every 100th variant
		if (!isset($target) && !$full_wgs && $count++ % 100 !== 0)
		{
			continue;
		}

		fwrite($handle_out, implode("\t", [ "chr{$chr}", $pos, $pos, strtoupper($ref), strtoupper($alt), $id ])."\n");
	}
	gzclose($handle);
	fclose($handle_out);
}

// annotate B-allele frequencies with BAM data
$annotated_variants = $parser->tempFile(".tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $filtered_variants -bam $in -out $annotated_variants -depth -name sample1", true);
if ($is_somatic)
{
	$annotated_variants_2 = $parser->tempFile(".tsv");
	$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $annotated_variants -bam $n_in -out $annotated_variants_2 -depth -name sample2", true);
	$annotated_variants = $annotated_variants_2;
}

// prepare IGV/SEG header
$seg_firstline = "#track graphtype=points viewLimits=-0.2:1.2 maxHeightPixels=80:80:80";
$seg_header = [ "Chromosome", "Start", "End", "Feature" ];
if ($is_somatic)
{
	$seg_header[] = "{$ps_name} BAFs (tumor)";
	$seg_header[] = "{$nor_name} BAFs (normal)";
}
else
{
	$seg_header[] = "{$ps_name} BAFs";
}
// write header lines
$handle_out = fopen($out, "w");
fwrite($handle_out, $seg_firstline . "\n");
fwrite($handle_out, implode("\t", $seg_header) . "\n");

// filter B-allele frequencies by depth
$handle = fopen($annotated_variants, "r");
while (!feof($handle))
{
	$row = nl_trim(fgets($handle));

	// skip empty lines, comments and header
	if ($row === "" || starts_with($row, "#"))
	{
		continue;
	}

	// skip low depth entries
	$cols = explode("\t", $row);
	if ($cols[7] < $min_dp || ($is_somatic && $cols[9] < $min_dp))
	{
		continue;
	}
	
	$row_out = [ $cols[0], $cols[1] - 1, $cols[2], $cols[5], $cols[6] ];
	if ($is_somatic)
	{
		$row_out[] = "${cols[8]}";
	}
	fwrite($handle_out, implode("\t", $row_out)."\n");
}

?>