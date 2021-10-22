<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("db_converter_dbscsnv", "Converts the dbscCNV flat files to VCF.GZ.");
// optional
$parser->addInfileArray("in",  "Input files as provided by dbscSNV.", false, false);
$parser->addInt("c_chr",  "Chromosome column index (0-based).", false);
$parser->addInt("c_pos",  "Chromosomal coordinate column index (0-based).", false);
$parser->addOutfile("out",  "Output file in VCF.GZ format.", false, false);
extract($parser->parse($argv));

// write VCF header
$tmp_file = temp_file(".vcf", "dbscsnv");
$out_fp = fopen2($tmp_file, "w");
fwrite($out_fp, "##fileformat=VCFv4.2\n");
fwrite($out_fp, "##fileDate=".date("Ymd")."\n");
fwrite($out_fp, "##INFO=<ID=ADA,Number=1,Type=Float,Description=\"dbscSNV ADA score.\">\n");
fwrite($out_fp, "##INFO=<ID=RF,Number=1,Type=Float,Description=\"dbscSNV RF score.\">\n");
fwrite($out_fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

foreach($in as $filename)
{
	$in_fp = fopen2($filename, "r");

	//parse input file 
	while(!feof($in_fp))
	{
		$line = trim(fgets($in_fp));
		if ($line=="" || starts_with($line, "chr\t")) continue;
		
		//parse line
		$parts = explode("\t", $line);
		if (count($parts)!=18)
		{
			trigger_error("Invalid number of columns (".count($parts).")!", E_USER_ERROR);
		}
		
		$chr = "chr".trim($parts[$c_chr]);
		$pos = trim($parts[$c_pos]);
		$ref = trim($parts[2]);
		$alt = trim($parts[3]);
		$ada_score = trim($parts[16]);
		if (is_numeric($ada_score)) $ada_score = number_format($ada_score, 3);
		$rf_score = trim($parts[17]);
		if (is_numeric($rf_score)) $rf_score = number_format($rf_score, 3);
		fwrite($out_fp, "$chr\t$pos\t.\t$ref\t$alt\t.\t.\tADA={$ada_score};RF={$rf_score}\n");
	}
	fclose($in_fp);
}
fclose($out_fp);


//Normalize file and sort
$parser->exec(get_path("ngs-bits")."/VcfSort", " -in {$tmp_file} -out {$out}");

?>