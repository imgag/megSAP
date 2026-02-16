<?php

/** 
	@page baf_germline
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("baf_germline", "Generate SEG file that contains b-allele frequencies for a germline sample.");
$parser->addString("name",  "Sample name used in output header.", false);
$parser->addInfile("vcf", "Variant list with SNPs to use (in VCF/VCF.GZ format).", false);
$parser->addOutfile("out",  "Output IGV file.", false);
//optional
$parser->addInt("min_dp", "Minimum depth of SNP locations in BAM.", true, 20);
$parser->addInt("downsample", "Enable downsampling, i.e. only every n-th SNP is used to calculate BAFs.", true, 0);
extract($parser->parse($argv));

//write output file header
$handle_out = fopen2($out, "w");
fwrite($handle_out, "#track graphtype=points viewLimits=-0.2:1.2 maxHeightPixels=80:80:80\n");
$seg_header = [ "Chromosome", "Start", "End", "Feature", "{$name} BAF" ];
fwrite($handle_out, implode("\t", $seg_header) . "\n");

//extract SNP list
$handle = gzopen2($vcf, "r");
$snps_passed = 0;
$snps_used = 0;
$done = [];
while (!feof($handle))
{
	$line = trim(fgets($handle));

	// skip empty lines and comments
	if ($line=="" || $line[0]=="#") continue;

	$cols = explode("\t", $line);
	if (count($cols)!=10) trigger_error("VCF file contains not exactly 10 columns:\n$line", E_USER_ERROR);

	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = $cols;

	// skip indels and multi-allelic sites
	if (strlen($ref) > 1 || strlen($alt) > 1) continue;
	
	// skip gene conversion events (Dragen)
	if (contains($info, 'EVENTTYPE=GENE_CONVERSION')) continue;
	
	++$snps_passed;
	
	// skip same positions
	$tag = $chr."_".$pos;
	if (isset($done[$tag])) continue;
	$done[$tag] = true;
	
	//extract DP from sample
	$sample = explode(":", $sample);
	$format = explode(":", $format);
	$sample = array_combine($format, $sample);
	$dp = "";
	if (isset($sample['DP']))
	{
		$tmp = $sample['DP'];
		if (is_numeric($tmp))
		{
			$dp = $tmp;
		}
	}
	if ($dp=="") continue;
	
	//apply depth cutoff
	if ($dp<$min_dp) continue;
	
	//extract AF from sample
	$af = "";
	if (isset($sample['AF']))
	{
		$tmp = $sample['AF'];
		if (is_numeric($tmp)) //skip faulty variants with multiple values caused by vcfallelicprimitves
		{
			$af = $tmp;
		}
	}
	if (isset($sample['AO'])) //freebayes does not report AF, so we need to calculate it from AO
	{
		$tmp = $sample['AO'];
		if (is_numeric($tmp)) //skip faulty variants with multiple values caused by vcfallelicprimitves
		{
			$af = $tmp/$dp;
		}
	}
	if (isset($sample['FREQ'])) //varscan reports AF as FREQ, e.g. '99.62%'
	{
		$tmp = strtr($sample['FREQ'], ['%'=>'']);
		if (is_numeric($tmp)) //skip faulty variants with multiple values caused by vcfallelicprimitves
		{
			$af = $tmp/100.0;
		}
	}
	if (isset($sample['VAF'])) //deepvariant/deepsomatic reports AF as VAF
	{
		$tmp = $sample['VAF'];
		if (is_numeric($tmp)) //skip faulty variants with multiple values caused by vcfallelicprimitves
		{
			$af = $tmp;
		}
	}
	if ($af=="") continue;
	
	// perform downsampling when needed
	if ($downsample>0 && $snps_passed%$downsample!=0) continue;
	++$snps_used;
	
	//write output
	$row_out = array($chr, $pos-1, $pos, ".", number_format($af, 4));
	fwrite($handle_out, implode("\t", $row_out)."\n");
}
gzclose($handle);
fclose($handle_out);

//print statistics
print "{$snps_passed} SNPs found in VCF.\n";
print "{$snps_used} SNPs used for BAF calculation.\n";

?>
