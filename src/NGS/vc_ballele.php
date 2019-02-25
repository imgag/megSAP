<?php

/**
	@page vc_ballele
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("vc_ballele", "Generate SEG file with b-allele frequencies from pre-computed variants.");
$parser->addInfile("vcf",  "Input VCF file.", false);
$parser->addInfile("bam",  "Input BAM file.", false);
$parser->addOutfile("out",  "Output IGV file.", false);

//optional
$parser->addInfile("n_bam",  "Input normal file in BAM format (somatic mode).", true);
$parser->addFlag("downsample", "Downsample input variant list (e.g. for WGS).");
$parser->addInt("min_dp", "Minimum depth of SNPs in BAM.", true, 20);
$parser->addInt("factor", "Factor used for downsampling.", true, 100);
$parser->addFlag("depth", "Add depth column(s) to output file.");
$parser->addString("build", "Genome build to use.", true, "GRCh37");
extract($parser->parse($argv));

$ps_name = basename($bam, ".bam");
$is_somatic = isset($n_bam);
if ($is_somatic)
{
	$nor_name = basename($n_bam, ".bam");
}

//downsample input variants if requested
if ($downsample)
{
	$downsampled_vcf = $parser->tempFile("_downsampled.vcf");
	$pipeline = [];
	$pipeline[] = [ "gzip", "-cd {$vcf}" ];
	$pipeline[] = [ "awk", "'/^#/ || NR%{$factor} == 0' > {$downsampled_vcf}" ]; //use awk to keep headers and every n-th line
	$parser->execPipeline($pipeline, "downsample vcf");
	$vcf = $downsampled_vcf;
}

//annotate variant frequencies from BAM
$annotated1 = $parser->tempFile("_sample1.tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $vcf -bam $bam -out $annotated1 -depth -name sample1 -ref ".genome_fasta($build), true);
if ($is_somatic)
{
	$annotated2 = $parser->tempFile("_sample2.tsv");
	$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $annotated1 -bam $n_bam -out $annotated2 -depth -name sample2 -ref ".genome_fasta($build), true);
}

//extract relevant columns and filter by depth
$tmp_baf = $parser->tempFile("baf.tsv");
$pipeline = [];
if ($is_somatic)
{
	$pipeline[] = [ get_path("ngs-bits")."/TsvSlice", "-in {$annotated2} -cols chr,start,end,ID,sample1_freq,sample2_freq,sample1_depth,sample2_depth" ];
	$pipeline[] = [ get_path("ngs-bits")."/TsvFilter", "-filter 'sample1_depth >= {$min_dp}'" ];
	$pipeline[] = [ get_path("ngs-bits")."/TsvFilter", "-filter 'sample2_depth >= {$min_dp}' -out {$tmp_baf}" ];
}
else
{
	$pipeline[] = [ get_path("ngs-bits")."/TsvSlice", "-in {$annotated1} -cols chr,start,end,ID,sample1_freq,sample1_depth" ];
	$pipeline[] = [ get_path("ngs-bits")."/TsvFilter", "-filter 'sample1_depth >= {$min_dp}' -out {$tmp_baf}" ];
}
//$pipeline[] = [ "awk", "'!/^#/' > {$tmp_baf}" ];
$parser->execPipeline($pipeline, "rewrite baf");

//write IGV/SEG header
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

$file_out_h = fopen($out, "w");
fwrite($file_out_h, $seg_firstline . "\n");
fwrite($file_out_h, implode("\t", $seg_header) . "\n");

//rewrite temporary BAF file, 0-based coordinates, drop depth columns if not specified
$file_in_h = fopen($tmp_baf, "r");
while(!feof($file_in_h)) {
	$line = nl_trim(fgets($file_in_h));
	if (starts_with($line, "#"))
	{
		continue;
	}
	$parts = explode("\t", $line);
	if (count($parts) >= 6)
	{
		//rewrite start coordinate
		$parts[1] = $parts[1] - 1;

		//remove depth column(s)
		if (!$depth)
		{
			array_pop($parts);
			if ($is_somatic)
			{
				array_pop($parts);
			}
		}

		fwrite($file_out_h, implode("\t", $parts) . "\n");
	}
}
fclose($file_in_h);
fclose($file_out_h);

?>
