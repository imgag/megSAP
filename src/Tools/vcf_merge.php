<?php
/** 
	@page vcf_merge
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("vcf_merge", "Merges single-sample VCF files to one multi-sample VCF file.");
$parser->addInfileArray("in", "Input VCF.GZ files.", false);
$parser->addOutfile("out", "Output VCF file.", false);
$parser->addInfile("roi", "Target region filter", true);
extract($parser->parse($argv));

//init
$ngsbits = get_path("ngs-bits");

//store variant sample associations
$samples = [];
$vars =  array();
foreach($in as $filename)
{
	print "Processing $filename...\n";
	flush();
	if (!ends_with($filename, ".vcf.gz"))
	{
		trigger_error("Input file $filename does have '.vcf.gz' extension.", E_USER_ERROR);
	}
	
	$command = "zcat {$filename}";
	if ($roi!="") $command .= " | {$ngsbits}/VcfFilter -reg {$roi}";

	$sample_current = "";
	list($file) = exec2($command);
	foreach($file as $line)
	{
		$line = trim($line);
		if($line=="") continue;
		
		//header
		if($line[0]=="#")
		{
			if (starts_with($line, "##")) continue;
			$sample_current = explode("\t", $line)[9];
			$samples[] = $sample_current;
			
			continue;
		}
		
		//content
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = explode("\t", $line);
		
		//check format column is same
		if ($format_exp=="") $format_exp = $format;
		if (!starts_with($format, "GT:")) trigger_error("Format column of {$filename} is not valid: '{$format}'", E_USER_ERROR);
		
		//store data
		$tag = "{$chr}\t{$pos}\t{$ref}\t{$alt}";
		$vars[$tag][$sample_current] = substr($sample, 0, 3);
	}
}

//write output
$h = fopen2($out, 'w');
fputs($h, "##fileformat=VCFv4.2\n");
fputs($h, "##INFO=<ID=HOM,Number=1,Type=String,Description=\"Number of samples with homzygous variant.\">\n");
fputs($h, "##INFO=<ID=HET,Number=1,Type=String,Description=\"Number of samples with heterozygous variant.\">\n");
fputs($h, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
fputs($h, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".implode("\t", $samples)."\n");
foreach($vars as $tag => $sample2geno)
{
	list($chr, $pos, $ref, $alt) = explode("\t", $tag);
	
	//count het/hom
	$het = 0;
	$hom = 0;
	foreach($sample2geno as $sample => $geno)
	{
		$geno = strtr($geno, ["|"=>"/", "."=>"0"]);
		if ($geno=="1/1") ++$hom;
		else if ($geno=="1/0" || $geno=="0/1") ++$het;
		else trigger_error("Genotype column of variant {$chr}:{$pos} {$ref}>{$alt} in {$filename} is not valid: '{$geno}'", E_USER_ERROR);
	}
	fputs($h, "{$chr}\t{$pos}\t.\t{$ref}\t{$alt}\t30\tPASS\tHET={$het};HOM={$hom}\tGT");
	foreach($samples as $sample)
	{
		if(isset($sample2geno[$sample]))
		{
			fputs($h, "\t".$sample2geno[$sample]);
		}
		else
		{
			fputs($h, "\t0/0");
		}
	}
	fputs($h, "\n");
}
fclose($h);

?>