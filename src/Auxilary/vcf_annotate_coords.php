<?php
/** 
	@page vcf_annotate_coords
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command input_line arguments
$parser = new ToolBase("vcf_annotate_coords", "Annotates the VCF with variant info.");
$parser->addInfile("in", "Input VCF file.", false);
$parser->addOutfile("out", "Output VCF file.", false);
$parser->addString("name", "Field name in INFO column.", true, "hg19_var");
extract($parser->parse($argv));

// read VCF
$input_lines = file($in, FILE_IGNORE_NEW_LINES);
$output_lines = [];

// parse file
foreach ($input_lines as $input_line) 
{
	if ($input_line[0]=="#")
	{
		if(starts_with($input_line, "#CHROM"))
		{
			$output_lines[] = "##INFO=<ID={$name},Number=1,Type=String,Description=\"Variant info (position and nucleotide change) of variant\">";
		}
		$output_lines[] = $input_line;
	}
	else
	{
		$split_line = explode("\t", $input_line);
		$var_info_entry = $name."=".implode("_", array($split_line[0], $split_line[1], $split_line[3], $split_line[4]));
		if($split_line[7] == "." || $split_line[7] == "")
		{
			$split_line[7] = $var_info_entry;
		}
		else
		{
			$split_line[7].= ";".$var_info_entry;
		}
		$output_lines[] = implode("\t", $split_line);
	}
}

file_put_contents($out, implode("\n", $output_lines));

?>