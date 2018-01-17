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
$parser->addInt("min_af", "Minimum allele frequency of SNPs to use (in 1000g data).", true, 0.01);
$parser->addInt("min_dp", "Minimum depth of SNPs in BAM.", true, 20);
$parser->addInfile("test", "Use this prefiltered list of common SNVs in vcf format. VCF must contain an AF field that denotes the SNV's population frequency.", true);
extract($parser->parse($argv));

//init
$ps_name = basename($in, ".bam");
$is_wgs = empty($target) && empty($test);
$is_somatic = !is_null($n_in);
$nor_name = null;
if($is_somatic) $nor_name = basename($n_in,".bam");

//filter SNPs by ROI
if(empty($test))
{
	$known_variants = get_path("data_folder")."/dbs/1000G/1000g_v5b.vcf.gz";
	if(!$is_wgs)
	{
		$known_variants = $parser->tempFile(".vcf");
		$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in ".get_path("data_folder")."/dbs/1000G/1000g_v5b.vcf.gz"." -reg $target -out $known_variants", true);
	}
}
else
{
	$known_variants = $test;
}

//filter known variants (SNPs with high population AF)
$count = 0;
$handle = gzopen($known_variants, "r");
if($handle===FALSE) trigger_error("Could not open file $known_variants for reading.",E_USER_ERROR);
$filtered_variants = $parser->tempFile(".tsv");
$handle_out = fopen($filtered_variants, "w");
if($handle_out===FALSE) trigger_error("Could not open file $filtered_variants for writing.",E_USER_ERROR);
fwrite($handle_out, "#chr\tstart\tend\tref\tobs\tid\n");
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
	
	//WGS: use only every 100th variant (otherwise this filtering will take too long)
	++$count;
	if($is_wgs && $count%100!=0) continue;
	
	fwrite($handle_out, "chr$chr\t$pos\t$pos\t".strtoupper($ref)."\t".strtoupper($alt)."\t$id\n");
}
gzclose($handle);
fclose($handle_out);

//annotate BAFs and depth
$annotated_variants = $parser->tempFile(".tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $filtered_variants -bam $in -out $annotated_variants -depth -name sample1", true);
if($is_somatic)
{
	$tmp_file = $parser->tempFile(".tsv");
	$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $annotated_variants -bam $n_in -out $tmp_file -depth -name sample2", true);
	$annotated_variants = $tmp_file;
}

//write output (filter by depth)
$handle = fopen($annotated_variants, "r");
if($handle===FALSE) trigger_error("Could not open file $annotated_variants for reading.",E_USER_ERROR);
$handle_out = fopen($out, "w");
if($handle_out===FALSE) trigger_error("Could not open file $out for writing.",E_USER_ERROR);
fwrite($handle_out, "#track graphtype=points viewLimits=-0.2:1.2 maxHeightPixels=80:80:80\n");
$header = "Chromosome	Start	End	Feature";
if(is_null($n_in))	$header .= "	$ps_name";
else	$header .= "	$ps_name (tumor)	$nor_name (normal)";
fwrite($handle_out, $header."\n");
while(!feof($handle))
{
	$row = nl_trim(fgets($handle));
	if ($row=="" || starts_with($row, "#")) continue;
	
	$cols = explode("\t", $row);
	if($cols[7]<$min_dp) continue;
	if($is_somatic && $cols[9]<$min_dp) continue;
	
	//#track graphtype=points viewLimits=-1.5:1.5 maxHeightPixels=80:80:80
	//Chromosome	Start	End	Feature	Tumor	Normal
	//chr1	2488152	2488153	SNP1	0.5	1
	//chr1	2491305	2491306	SNP2	-0.5	1
	//chr1	2494329	2494330	SNP3	-0.1267	1
	//chr1	4849383	4849384	SNP4	-0.4938	1
	
	$row_out = "${cols[0]}\t".($cols[1]-1)."\t${cols[2]}\t${cols[5]}\t${cols[6]}";
	if($is_somatic) 
	{
		$row_out .= "\t${cols[8]}";
	}
	fwrite($handle_out, $row_out."\n");
}
fclose($handle);
fclose($handle_out);

?>
