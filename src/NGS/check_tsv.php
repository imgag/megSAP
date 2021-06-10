<?php

/** 
	@page check_tsv
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("check_tsv", "Performs TSV file format check.");
$parser->addInfile("in",  "Input file in TSV format.", false);
//optional
$parser->addInt("limit", "The number of variants to check for format-specific checks. '0' means all.", true, 1000);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//define format in the following format: column index, name, type
$col_desc = array(
	array(0, "chr", "chr"),
	array(1, "start", "num"),
	array(2, "end", "num"),
	array(3, "ref", "nuc"),
	array(4, "obs", "nuc"),
	//array(5, "genotype", "genotype")
	);

//load input
$i = -1;
$handle = fopen2($in, "r");
while(!feof($handle))
{
	++$i;
	$line = nl_trim(fgets($handle));
	
	//skip empty lines
	if($line=="") continue;
	
	//skip comment lines
	if (starts_with($line, "##")) continue;
	
	//check headers
	if (starts_with($line, "#"))
	{
		$headers = explode("\t", substr($line, 1));
		foreach ($col_desc as $col)
		{
			list ($col_index, $col_name) = $col;
			if ($headers[$col_index]!=$col_name)
			{
				trigger_error("Header with index '$col_index' should be '$col_name', but is '".$headers[$col_index]."'!", E_USER_ERROR);
			}
		}
		continue;
	}
	
	//content lines
	$row_values = explode("\t", $line);
	
	//check that number of columns is correct
	if (count($row_values)!=count($headers))
	{
		trigger_error("Column count is ".count($row_values).", but should be ".count($headers)." in line:\n$line", E_USER_ERROR);
	}
	
	//check column types
	foreach ($col_desc as $col)
	{
		list ($col_index, $col_name, $col_type) = $col;
		$value = $row_values[$col_index];
		
		switch ($col_type)
		{
			case "num":
				if($value!="" && $value!="." && !is_numeric($value))	
				{
					trigger_error("Invalid value '$value' for numeric column (index $col_index, name $col_name) in line ".($i).".", E_USER_ERROR);
				}
				break;
			case "chr":
				break;
			case "nuc":
				if(	$value!="-" && !preg_match("/^[ACTGN]+$/i", $value)) 
				{
					trigger_error("Invalid value '$value' for nucleic acid column (index $col_index, name $col_name) in line ".($i).".", E_USER_ERROR);
				}
				break;
			case "genotype":
				if($value!="hom" && $value!="het" && $value!=".")
				{
					trigger_error("Invalid value '$value' for genotype column (index $col_index, name $col_name) in line ".($i).".", E_USER_ERROR);
				}
				break;
			default:
				trigger_error("Unknown column type '$col_type'!", E_USER_ERROR);
		}
	}
	
	//check TSV format
	if ($limit<=0 || $i<=$limit)
	{
		list($chr, $start, $end, $ref, $obs) =  $row_values;

		//check SNP
		$is_snp = $ref!="-" && $obs!="-" && strlen($ref)==1 && strlen($obs)==1;
		if($is_snp)
		{
			//one-based position format
			if($start!=$end)
			{
				trigger_error("SNP position end format is not one-based in line with index '$i'.", E_USER_ERROR);
			}
				
			//check reference genome
			$ref_seq = get_ref_seq($build, $chr, $start, $end);
			if(strcasecmp($ref_seq, $ref)!=0)
			{
				trigger_error("SNP ref base '$ref' does not match genome sequence '$ref_seq' in line with index '$i'.", E_USER_ERROR);
			}
		}
		//check deletions
		else if($ref!="-")
		{
			//check reference genome
			$ref_seq = get_ref_seq($build, $chr, $start, $end);
			if(strcasecmp($ref_seq, $ref)!=0)
			{
				trigger_error("INDEL reference sequence '$ref' does not match genome sequence '$ref_seq' in line with index '$i'.", E_USER_ERROR);
			}	
		}
	}
}

fclose($handle);

?>
