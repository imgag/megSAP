<?php

/** 
	@page baf_germline
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("baf_germline", "Generate SEG file that contains b-allele frequencies for a germline sample.");
$parser->addInfile("bam",  "Input BAM file.", false);
$parser->addInfile("vcf", "Variant list with SNPs to use (in VCF/VCF.GZ format).", false);
$parser->addOutfile("out",  "Output IGV file.", false);
//optional
$parser->addInt("min_dp", "Minimum depth of SNP locations in BAM.", true, 20);
$parser->addFlag("depth", "Add depth column(s) to 'out'.", true, 20);
$parser->addInt("downsample", "Enable downsampling, i.e. only every n-th SNP is used to calculate BAFs.", true, 0);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
extract($parser->parse($argv));

$ps_name = basename($bam, ".bam");

//extract SNP list
$snps_filtered = $parser->tempFile("snps.tsv");

$handle = gzopen2($vcf, "r");

$handle_out = fopen2($snps_filtered, "w");
fwrite($handle_out, "#chr\tstart\tend\tref\tobs\n");

$snps_passed = 0;
$snps_used = 0;
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

	// skip indels and multi-allelic sites
	if (strlen($ref) > 1 || strlen($alt) > 1)
	{
		continue;
	}
	
	++$snps_passed;
	
	if ($downsample>0 && $snps_passed%$downsample!=0)
	{
		continue;
	}

	++$snps_used;	
	fwrite($handle_out, implode("\t", [ $chr, $pos, $pos, strtoupper($ref), strtoupper($alt) ])."\n");
}
gzclose($handle);
fclose($handle_out);
print "{$snps_passed} SNPs found in VCF.\n";
print "{$snps_used} SNPs used for BAF calculation.\n";

//annotate B-allele frequencies from BAM
$annotated_variants = $parser->tempFile(".tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $snps_filtered -bam $bam -out $annotated_variants -depth -name sample1 -ref ".genome_fasta($build), true);

//prepare IGV/SEG header
$seg_firstline = "#track graphtype=points viewLimits=-0.2:1.2 maxHeightPixels=80:80:80";
$seg_header = [ "Chromosome", "Start", "End", "Feature", "{$ps_name} BAF" ];
if ($depth) $seg_header[] = "{$ps_name} DP";

//write header lines
$non_unique_out = $parser->tempFile(".vcf");
$handle_out = fopen2($non_unique_out, "w");
fwrite($handle_out, $seg_firstline . "\n");
fwrite($handle_out, implode("\t", $seg_header) . "\n");

//filter B-allele frequencies by depth
$handle = fopen2($annotated_variants, "r");
while (!feof($handle))
{
	$row = nl_trim(fgets($handle));

	// skip empty lines, comments and header
	if (trim($row)=="" || $row[0]=="#")
	{
		continue;
	}
	
	$cols = explode("\t", $row);
	list($chr, $start, $end, $ref, $obs, $af, $dp) = $cols;

	// skip low depth entries
	if ($dp<$min_dp) continue;
	
	$row_out = array($chr, $start-1, $end, ".", $af);
	if ($depth)
	{
		$row_out[] = $dp;
	}
	fwrite($handle_out, implode("\t", $row_out)."\n");
}
fclose($handle_out);

//Make sure, out file only contains unique variants
$parser->exec("uniq","$non_unique_out $out",true);

?>
