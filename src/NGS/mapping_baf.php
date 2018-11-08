<?php

/** 
	@page mapping_baf
	@todo use sample VCF as input to make the BAF file more informative - remove gnomAD dependency - move to 'vc' step.
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("mapping_baf", "Generate SEG file that contains b-allele frequencies from common snps.");
$parser->addInfile("in",  "Input BAM file.", false);
$parser->addOutfile("out",  "Output IGV file.", false);
//optional
$parser->addInfile("target",  "Target region BED file defining which population SNPs are extracted from 'snp_db'. If unspecified, WGS mode (downsampling to every 100th SNP) will be used.", true);
$parser->addInfile("n_in",  "Input normal file in BAM format (somatic mode).", true);
$parser->addInfile("sites", "Additional sites to use (VCF/VCF.GZ format).", true);
$parser->addFloat("min_af", "Minimum allele frequency of SNPs in 'snp_db' to use.", true, 0.01);
$parser->addInt("min_dp", "Minimum depth of SNPs in BAM.", true, 20);
$parser->addString("snp_db", "SNP database in VCF.GZ format from which population SNPs are extracted. Records must contain AF field.", true, get_path("data_folder")."/dbs/gnomAD/gnomAD_genome_r2.0.2.vcf.gz");
$parser->addFlag("depth", "Add depth column(s) to 'out'.", true, 20);
extract($parser->parse($argv));

//check if base VCF exists
if (!file_exists($snp_db))
{
	trigger_error("SNP database '{$snp_db}' does not exist!", E_USER_ERROR);
}

$ps_name = basename($in, ".bam");
$is_somatic = isset($n_in);
if ($is_somatic)
{
	$nor_name = basename($n_in, ".bam");
}

// if WGS flag is set, target region is not relevant
$is_wgs = !isset($target);
if ($is_wgs)
{
	$shortname = "WGS";
}
else
{
	$shortname = basename($target, ".bed") . "_" . sha1_file($target);
}

//use pre-computed SNP list if possible
$filtered_variants = get_path("local_data")."/baf_".basename($snp_db, ".vcf.gz")."_af{$min_af}_{$shortname}.tsv";
if (!file_exists($filtered_variants))
{
	trigger_error("Filtered variant list does not exist, creating at '{$filtered_variants}'.", E_USER_NOTICE);

	// filter base VCF by target region (only for enrichment)
	if ($is_wgs)
	{
		$tmp_filtered_by_region = $snp_db;
	}
	else
	{
		$tmp_filtered_by_region = $parser->tempFile("_filtered_by_region.vcf");
		$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in {$snp_db} -reg {$target} -out {$tmp_filtered_by_region}", true);
	}

	// filter known variants (SNPs with high population AF)
	$count = 0;

	$handle = gzopen($tmp_filtered_by_region, "r");
	if ($handle === FALSE)
	{
		trigger_error("Could not open file '{$snp_db}' for reading.", E_USER_ERROR);
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

		// skip indels and multi-allelic sites
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
		if ($is_wgs && $count++ % 100 !== 0)
		{
			continue;
		}

		fwrite($handle_out, implode("\t", [ "chr{$chr}", $pos, $pos, strtoupper($ref), strtoupper($alt), $id ])."\n");
	}
	gzclose($handle);
	fclose($handle_out);
}

//add sites from provided vcf
if (isset($sites))
{
	$filtered_variants_add = $parser->tempFile("add.tsv");
	$parser->copyFile($filtered_variants, $filtered_variants_add);

	$handle = gzopen($sites, "r");
	$handle_out = fopen($filtered_variants_add, "a");
	while (!feof($handle))
	{
		$line = nl_trim(fgets($handle));
		//skip empty lines and comments
		if ($line == "" || starts_with($line, "#"))
		{
			continue;
		}

		$cols = explode("\t", $line, 8);
		if (count($cols) < 8)
		{
			trigger_error("VCF file line contains less than 8 columns:\n$line", E_USER_ERROR);
		}

		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = $cols;
		if (!($filter === "." || $filter === "PASS") ||	//skip filtered variants
			!chr_check($chr, 22, false) ||				//skip special chromosomes
			strlen($ref) > 1 || strlen($alt) > 1)		//skip indels, multi-allelic sites
		{
			continue;
		}

		fwrite($handle_out, implode("\t", [ $chr, $pos, $pos, strtoupper($ref), strtoupper($alt), $id ])."\n");
	}
	fclose($handle_out);

	//sort variants from SNP DB + added sites
	$filtered_variants_add_sorted = $parser->tempFile("sorted.tsv");
	$parser->exec(get_path("ngs-bits")."/VcfSort", "-in $filtered_variants_add -out $filtered_variants_add_sorted", true);

	$filtered_variants = $filtered_variants_add_sorted;
}

// annotate B-allele frequencies from BAM
$annotated_variants = $parser->tempFile(".tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $filtered_variants -bam $in -out $annotated_variants -depth -name sample1 -ref ".get_path("local_data")."/GRCh37.fa", true);
if ($is_somatic)
{
	$annotated_variants_2 = $parser->tempFile(".tsv");
	$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $annotated_variants -bam $n_in -out $annotated_variants_2 -depth -name sample2 -ref ".get_path("local_data")."/GRCh37.fa", true);
	$annotated_variants = $annotated_variants_2;
}

// prepare IGV/SEG header
$seg_firstline = "#track graphtype=points viewLimits=-0.2:1.2 maxHeightPixels=80:80:80";
$seg_header = [ "Chromosome", "Start", "End", "Feature" ];
if ($is_somatic)
{
	$seg_header[] = "{$ps_name} BAF (tumor)";
	$seg_header[] = "{$nor_name} BAF (normal)";
}
else
{
	$seg_header[] = "{$ps_name} BAF";
}
if ($depth)
{
	if ($is_somatic)
	{
		$seg_header[] = "{$ps_name} DP (tumor)";
		$seg_header[] = "{$nor_name} DP (normal)";
	}
	else
	{
		$seg_header[] = "{$ps_name} DP";
	}
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
	if (trim($row)=="" || $row[0]=="#")
	{
		continue;
	}
	
	$cols = explode("\t", $row);
	list($chr, $start, $end, $ref, $obs, $id, $af, $dp) = $cols;

	// skip low depth entries
	if ($is_somatic)
	{
		$af_nor = $cols[8];
		$dp_nor = $cols[9];
	}
	if ($dp<$min_dp || ($is_somatic && $dp_nor < $min_dp))
	{
		continue;
	}
	
	$row_out = array($chr, $start-1, $end, $id, $af);
	if ($is_somatic) $row_out[] = $af_nor;
	if ($depth)
	{
		$row_out[] = $dp;
		if ($is_somatic) $row_out[] = $dp_nor;
	}
	fwrite($handle_out, implode("\t", $row_out)."\n");
}
fclose($handle_out);

?>
