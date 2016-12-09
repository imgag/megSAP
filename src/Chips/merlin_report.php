<?php

/** 
	@page merlin_report
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");


// parse command line arguments
$parser = new ToolBase("merlin_report", "\$Rev: 789 $", "Converts linkage region from a merlin parametric model analysis to a BED file.");
$parser->addInfile("in",  "Input MODEL file", false);
$parser->addOutfile("out",  "Output BED file", true);
extract($parser->parse($argv));

//retruns the dbSNP identifier where linkage starts or ends from a merlin model file
function extract_rs_numbers($merlin_tbl_file_name,  $cutoff)
{
	$output = array();
	
	$linkage=FALSE;
	$last_entry = "";	
	$old_chrom = "";
	
	$merlin_tbl_file = file($merlin_tbl_file_name);
	array_shift($merlin_tbl_file); //skip header
	foreach($merlin_tbl_file as $line) 
	{ 
		$parts = explode("\t", $line);
		if (count($parts)<5) continue;
		
		list($current_chrom, , $rsname, , $lod_score) = $parts;
		
		//check if chromsome ended with an open linkage block 
		if ($current_chrom!=$old_chrom)
		{
			if ($linkage)
			{
				//the linkage blocks ends with the end of the chromosome
				$output[] = $last_entry;
				$linkage = FALSE;
			}
		}
		
		//check if a new linkage blocks started
		if ($lod_score>=$cutoff)
		{
			if (!$linkage)
			{
				//start of a linkage block
				$output[] = $rsname;
				$linkage = TRUE;
			}
		}
		
		//check if a linkage blocks ended
		else
		{
			if ($linkage)
			{
				//end of a linkage block should NOT be the first marker with low lod score, but the last with a high one
				$output[] = $last_entry;
				$linkage = FALSE;
			}
		}
		$last_entry = $rsname;
		$old_chrom = $current_chrom;
	}
	
	return array_map("trim", $output);
}

//Returns a dictionary rs-number => array(chr, pos)
function extract_dbsnp_entries($rs_numbers)
{
	//find dbSNP entries based on rsnumbers, includig those which only partial match
	$dbsnp_hits = array();
	exec("zegrep \"".implode("|", $rs_numbers)."\" ".get_path("data_folder")."/dbs/dbSNP/dbsnp_b147.vcf.gz", $dbsnp_hits);
		
	//check which grep results match the rs-numbers exactly and save these
	$output = array();
	foreach($dbsnp_hits as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		
		list($chr, $pos, $rs) = explode("\t", $line);		
		if (in_array($rs, $rs_numbers))
		{
			$output[$rs] = array("chr".$chr, $pos);
		}
	}
	
	return $output;
}

//build bedfile of linkage blocks, based on a file with rsnumbers and a file with their corresponding dbSNP entries, including (unwanted) partial matches
function build_bedfile($rs_numbers, $dbsnp_entries, $output_filename)
{
	//write bed file
	$line = "";
	$output = array();
	foreach ($rs_numbers as $rsnum)
	{		
		if ($line=="") //start of block (1. line)
		{
			$line = $dbsnp_entries[$rsnum][0]."\t".$dbsnp_entries[$rsnum][1];
		}
		else //end of block (2. line)
		{
			$output[] = $line."\t".$dbsnp_entries[$rsnum][1];
			$line = "";
		}
	}
	file_put_contents($output_filename, implode ("\n", $output));
}

$rs_numbers = extract_rs_numbers($in, 1.0);
$dbsnp_entries = extract_dbsnp_entries($rs_numbers);
build_bedfile($rs_numbers, $dbsnp_entries, $out);

?>