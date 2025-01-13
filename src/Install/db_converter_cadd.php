<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("db_converter_cadd", "Converts the CADD flat files (tsv.gz) to VCF(.gz).");
// optional
$parser->addInfile("in",  "Input file in tsv-format with CADD scores. ('-' for STDIN)", false, false);
$parser->addOutfile("out",  "Output file in vcf.gz-format with CADD scores. ('-' for STDOUT)", false, false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

// open input file
if ($in == "-") $in = "php://stdin";
$in_fp = fopen2($in, "r");

// open output file
if ($out == "-") $out = "php://stdout";
$out_fp = fopen2($out, "w");

// write VCF header
fwrite($out_fp, "##fileformat=VCFv4.2\n");
fwrite($out_fp, "##fileDate=".date("Ymd")."\n");
fwrite($out_fp, "##source={$in}\n");
fwrite($out_fp, "##reference=".genome_fasta($build, false)."\n");
fwrite($out_fp, "##INFO=<ID=CADD,Number=1,Type=Float,Description=\"CADD PHRED score of this variant.\">\n");
fwrite($out_fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

// parse input file 
while(!feof($in_fp))
{
	$line = fgets($in_fp);
	// skip comments/header
	if (starts_with($line, "#")) continue;
	// skip empty lines
	if (trim($line) == "") continue;
	
	// parse line
	$split_line = explode("\t", $line);
	if (count($split_line) != 6)
	{
		trigger_error("Invalid number of columns (".count($split_line).") in line: $line", E_USER_ERROR);
	}
	$chr = "chr".chr_trim($split_line[0]);
	$pos = (int) $split_line[1];
	$ref = trim($split_line[2]);
	$alt = trim($split_line[3]);
	$cadd_score = trim($split_line[5]);

	// correct InDels in which the base after the InDel is given
	if (($ref[0] != $alt[0]) && (substr($ref, -1) == substr($alt, -1)) && ($pos > 1))
	{
		// correct position
		$pos -= 1;

		// get previous base
		$prev_base = get_ref_seq($build, $chr, $pos, $pos, 100000, false);

		// remove last base and add base in front:
		$ref = $prev_base.substr($ref, 0, -1);
		$alt = $prev_base.substr($alt, 0, -1);	
	}

	// write line
	fwrite($out_fp, "$chr\t$pos\t.\t$ref\t$alt\t.\t.\tCADD=$cadd_score\n");
}

fclose($in_fp);
fclose($out_fp);
