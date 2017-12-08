<?php

/** 
	@page mapping_bwa
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_baf", "Generate seg file that contains bafs.");
$parser->addInfile("bam",  "Input file in bam format.", false);
$parser->addInfile("target",  "Input file in bed format.", true);
$parser->addOutfile("out",  "Output file in seg format.", false);
//optional
extract($parser->parse($argv));

$sid = basename($bam,".bam");

// filter 1000g for target region and population freqeuency
$snps =  get_path("data_folder")."/dbs/1000G/1000g_v5b.vcf.gz";
$tmp_variants1 = $parser->tempFile("_variants_filtered.vcf");
$parser->exec(get_path("ngs-bits")."/VariantFilterRegions", "-in $snps -reg $target -out $tmp_variants1", true);

$min_af = 0.01;
$tmp_variants2 = $parser->tempFile("_variants_filtered.tsv");
$handle = fopen($tmp_variants1, "r");
$handle_out = fopen($tmp_variants2, "w");
$in_header = true;
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="" || starts_with($line, "#")) continue;
	
	$cols = explode("\t", $line);
	if (count($cols)<10) trigger_error("VCF file line contains less than 10 columns:\n$line", E_USER_ERROR);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = $cols;
	
	$skip_variant = true;
	$ii = explode(";",$info);
	foreach($ii as $i)
	{
		if(!contains($i,"="))	continue;
		list($n,$v) = explode("=",$i);
		if($n=="AF" && $v>=$min_af)	$skip_variant = false;
	}
	if($skip_variant)	continue;
	
	//parse data from VCF
	if(chr_check($chr, 22, false) === FALSE) continue; //skip bad chromosomes
	if(strlen($ref)>1 || strlen($alt)>1)	continue;	//skip indels
	
	$start = $pos;
	$end = $pos;
	$ref = strtoupper($ref);
	$alt = strtoupper($alt);
	fwrite($handle_out, "chr$chr\t$start\t$end\t$ref\t$alt\n");
}
fclose($handle);
fclose($handle_out);

$tmp_variants3 = $parser->tempFile("_annotated.tsv");
$parser->exec(get_path("ngs-bits")."/VariantAnnotateFrequency", "-in $tmp_variants2 -bam $bam -out $tmp_variants3 -depth", true);

$min_depth = 20;
$tmp_variants4 = $parser->tempFile("_variants_".$min_depth."x.vcf");
$handle = fopen($tmp_variants3, "r");
$handle_out = fopen($tmp_variants4, "w");
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="" || starts_with($line, "#")) continue;
	
	$cols = explode("\t", $line);
	if($cols[6]<$min_depth)	continue;
	
	fwrite($handle_out, implode("\t",$cols)."\n");
}
fclose($handle);
fclose($handle_out);

$handle = fopen($tmp_variants4, "r");
$handle_out = fopen($out, "w");
fwrite($handle_out, "#type=GENE_EXPRESSION\n");
fwrite($handle_out, "#track graphtype=points name=\"".$sid." CN baf\" midRange=0.25:0.75 color=0,0,255 altColor=255,0,0 viewLimits=0:1 maxHeightPixels=80:80:80\n");
fwrite($handle_out, "ID	chr	start	end	z-score\n");
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	
	if ($line=="") continue;
	if (starts_with($line, "#")) 	continue;

	$cols = explode("\t", $line);	
	fwrite($handle_out, $sid."\t${cols[0]}\t${cols[1]}\t${cols[2]}\t${cols[5]}\n");
}
fclose($handle);
fclose($handle_out);

?>
