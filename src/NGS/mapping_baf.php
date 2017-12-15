<?php

/** 
	@page mapping_baf
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("mapping_baf", "Generate SEG file that contains bafs.");
$parser->addInfile("in",  "Input BAM file.", false);
$parser->addOutfile("out",  "Output SEG file.", false);
//optional
$parser->addInfile("system",  "Processing system INI file (determined from 'in' by default).", true);
$parser->addInt("min_af", "Minimum allele frequency of SNPs to use (in 1000g data).", true, 0.01);
$parser->addInt("min_dp", "Minimum depth of SNPs in BAM.", true, 20);
extract($parser->parse($argv));

//init
$ps_name = basename($in, ".bam");
$sys = load_system($system, $ps_name);
$is_wgs = $sys['type']=="WGS";

//filter SNPs by ROI
if ($is_wgs)
{
	$known_variants = get_path("data_folder")."/dbs/1000G/1000g_v5b.vcf.gz";
}
else if($sys['target_file']!="")
{
	$known_variants = $parser->tempFile(".vcf");
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in ".get_path("data_folder")."/dbs/1000G/1000g_v5b.vcf.gz"." -reg ".$sys['target_file']." -out $known_variants", true);
}
else
{
	trigger_error("Calclating BAFs requires WGS sample or target region!", E_USER_ERROR); 
}

//filter known variants (SNPs with high AF)
$count = 0;
$handle = gzopen($known_variants, "r");
$filtered_variants = $parser->tempFile(".tsv");
$handle_out = fopen($filtered_variants, "w");
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="" || starts_with($line, "#")) continue;
	
	$cols = explode("\t", $line);
	if (count($cols)<8) trigger_error("VCF file line contains less than 8 columns:\n$line", E_USER_ERROR);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = $cols;
	
	//skip special chromosomes
	if(chr_check($chr, 22, false) === FALSE) continue;
	
	//skip indels
	if(strlen($ref)>1 || strlen($alt)>1)	continue;
	
	//skip low AF variants
	$af = 0.0;
	$info = explode(";", $info);
	foreach($info as $entry)
	{
		if (starts_with($entry, "AF="))
		{
			$af = substr($entry, 3);
		}
	}
	if ($af<$min_af) continue;
	
	//WGS: use only every 100th variant (otherwise this takes too long)
	++$count;
	if($is_wgs && $count%100!=0) continue;

	fwrite($handle_out, "chr$chr\t$pos\t$pos\t".strtoupper($ref)."\t".strtoupper($alt)."\n");
}
gzclose($handle);
fclose($handle_out);

//annotate BAFs and depth
$annotated_variants = $parser->tempFile(".tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $filtered_variants -bam $in -out $annotated_variants -depth", true);

//write output (filter by depth)
$handle_out = fopen($out, "w");
fwrite($handle_out, "#type=GENE_EXPRESSION\n");
fwrite($handle_out, "#track graphtype=points name=\"{$ps_name} BAF\" midRange=0.25:0.75 color=0,0,255 altColor=255,0,0 viewLimits=0:1 maxHeightPixels=80:80:80\n");
fwrite($handle_out, "ID	chr	start	end	z-score\n");
$handle = fopen($annotated_variants, "r");
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="" || starts_with($line, "#")) continue;
	
	$cols = explode("\t", $line);
	if($cols[6]<$min_dp) continue;
	
	fwrite($handle_out, $ps_name."\t${cols[0]}\t${cols[1]}\t${cols[2]}\t${cols[5]}\n");
}
fclose($handle);
fclose($handle_out);

?>
