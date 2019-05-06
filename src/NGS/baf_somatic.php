<?php

/** 
	@page baf_somatic
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("baf_somatic", "Generate SEG file that contains b-allele frequencies for a tumor-normal pair.");
$parser->addInfile("bam_t",  "Tumor BAM file.", false);
$parser->addInfile("bam_n",  "Normal BAM format.", false);
$parser->addInfile("vcf", "Variant list with SNPs to use (in VCF/VCF.GZ format).", false);
$parser->addOutfile("out",  "Output IGV file.", false);
//optional
$parser->addInt("min_dp", "Minimum depth of SNP locations in BAMs.", true, 20);
$parser->addFlag("depth", "Add depth column(s) to 'out'.", true, 20);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
extract($parser->parse($argv));


$name_t = basename($bam_t, ".bam");
$name_n = basename($bam_n, ".bam");

//extract SNP list
$snps_filtered = $parser->tempFile("snps.tsv");

$handle = gzopen($vcf, "r");
if ($handle === FALSE) trigger_error("Could not open file '{$vcf}' for reading.", E_USER_ERROR);

$handle_out = fopen($snps_filtered, "w");
if ($handle_out === FALSE) trigger_error("Could not open file '{$snps_filtered}' for writing.", E_USER_ERROR);
fwrite($handle_out, "#chr\tstart\tend\tref\tobs\n");

$snps_passed = 0;
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
	
	fwrite($handle_out, implode("\t", [ $chr, $pos, $pos, strtoupper($ref), strtoupper($alt) ])."\n");
}
gzclose($handle);
fclose($handle_out);
print "{$snps_passed} SNPs found in VCF.\n";

//annotate B-allele frequencies from BAM
$annotated_variants = $parser->tempFile(".tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $snps_filtered -bam $bam_t -out $annotated_variants -depth -name sample1 -ref ".genome_fasta($build), true);
$annotated_variants_2 = $parser->tempFile(".tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $annotated_variants -bam $bam_n -out $annotated_variants_2 -depth -name sample2 -ref ".genome_fasta($build), true);

//prepare IGV/SEG header
$seg_firstline = "#track graphtype=points viewLimits=-0.2:1.2 maxHeightPixels=80:80:80";
$seg_header = [ "Chromosome", "Start", "End", "Feature", "{$name_t} BAF (tumor)", "{$name_n} BAF (normal)" ];
if ($depth)
{
	$seg_header[] = "{$name_t} DP (tumor)";
	$seg_header[] = "{$name_n} DP (normal)";
}

//write header lines
$non_unique_out = $parser->tempFile(".vcf");
$handle_out = fopen($non_unique_out, "w");
fwrite($handle_out, $seg_firstline . "\n");
fwrite($handle_out, implode("\t", $seg_header) . "\n");

//filter B-allele frequencies by depth
$handle = fopen($annotated_variants_2, "r");
while (!feof($handle))
{
	$row = nl_trim(fgets($handle));

	// skip empty lines, comments and header
	if (trim($row)=="" || $row[0]=="#")
	{
		continue;
	}
	
	$cols = explode("\t", $row);
	list($chr, $start, $end, $ref, $obs, $af_t, $dp_t, $af_n, $dp_n) = $cols;

	// skip low depth entries
	if ($dp_t<$min_dp || $dp_n < $min_dp) continue;
	
	$row_out = array($chr, $start-1, $end, ".", $af_t, $af_n);
	if ($depth)
	{
		$row_out[] = $dp_t;
		$row_out[] = $dp_n;
	}
	fwrite($handle_out, implode("\t", $row_out)."\n");
}
fclose($handle_out);

//Make sure, out file only contains unique variants
$parser->exec("uniq","$non_unique_out $out",true);

?>