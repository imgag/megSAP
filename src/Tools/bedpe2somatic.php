<?php
/** 
	@page bedpe2somatic
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("bedpe2somatic", "Expands data from INFO/FORMAT columns into single columns.");
$parser->addInfile("in", "Input BEDPE file.", false);
$parser->addString("tid", "Processed sample ID of tumor.", false);
$parser->addString("nid", "Processed sample ID of normal control.", false);
$parser->addOutfile("out", "Output file.", false);
extract($parser->parse($argv));


################################## functions ##################################
///Returns information from FORMAT column and corresponding data column as associative array
function expand_format_col($format_desc_col, $data_col)
{
	if(count($format_desc_col) != count($data_col))
	{
		trigger_error("Row count mismatch.", E_USER_ERROR);
	}
	
	$expanded_cols = array();

	for($row=0; $row<count($format_desc_col); ++$row)
	{
		$f_parts = explode(":",$format_desc_col[$row]); //parts of format description
		$d_parts = explode(":",$data_col[$row]);  //parts of corresponding data
		
		if(count($f_parts) != count($d_parts))
		{
			trigger_error("FORMAT count and col count differs.", E_USER_ERROR);
		}
		

		for($i=0;$i<count($f_parts); ++$i)
		{
			$expanded_cols[$f_parts[$i]][] = $d_parts[$i];
		}
	}
	
	$col_pr_ref = array();
	$col_pr_alt = array();
	$col_sr_ref = array();
	$col_sr_alt = array();

	for($row=0; $row<count($expanded_cols["PR"]); ++$row)
	{
		$pr_parts = explode(",", $expanded_cols["PR"][$row]);
		(count($pr_parts) == 2 )? list($pr_ref, $pr_alt) = $pr_parts : list($pr_ref, $pr_alt) = array(".", ".");
		
		$sr_parts = explode(",", $expanded_cols["SR"][$row]);
		(count($sr_parts) == 2)? list($sr_ref, $sr_alt) = $sr_parts : list($sr_ref, $sr_alt) = array(".", ".");
		
		$col_pr_ref[] = $pr_ref;
		$col_pr_alt[] = $pr_alt;
		
		$col_sr_ref[] = $sr_ref;
		$col_sr_alt[] = $sr_alt;
	}
	
	return array("PR_REF" => $col_pr_ref, "PR_ALT" => $col_pr_alt, "SR_REF" => $col_sr_ref, "SR_ALT" => $col_sr_alt);
}

///Returns description of info column as associative array 
function get_info_fields($comments)
{
	$info_col_description = array();
	foreach($comments as $comment)
	{
		if(!starts_with($comment, "#INFO=")) continue;
		
		$raw = trim($comment);
		//Remove initial "INFO=<" and ">" at the end of the string
		$raw =  substr($raw,7,strlen($raw)-8);
		
		$temp_data = array();
		
		//Special treatment for Description="" part, because it can contain additional ",", which is used as a separator
		$desc = substr($raw,strpos($raw,"Description="));
		$desc = substr($desc,12); //Remove "Description="
		$desc = substr($desc,1,strlen($desc)-2); //Remove '"' at beginning and end of Description 
		$temp_data["Description"] = $desc;
		
		//Parse remaining string without parts from "Description="-field
		$raw = substr($raw,0,strpos($raw,"Description=")-1);
		
		$parts = explode(",",$raw);
		
		$key = "";
		
		foreach($parts as $part)
		{
			list($id,$val) = explode("=",$part);
			if($id == "ID") 
			{
				$key = $val;
				continue;
			}
			$temp_data[$id] = $val;
		}
		
		$info_col_description[$key] = $temp_data;
	}

	return $info_col_description;
}

//Extracts certain fields from INFO column and returns results as associative array containing the colums.
//The function returns the keys "flags", "somaticscore" and "junction_somaticscore"
function expand_info_col($info_col)
{
	//Create list of all keys that occur in INFO column
	$all_keys = array();
	$row_count = 0;
	foreach($info_col as $field)
	{
		$parts = explode(";", $field);
		foreach($parts as $raw_part)
		{
			$pair = explode("=", $raw_part);
			
			$key = "";
			$value = "";
			
			if(count($pair) == 2)
			{
				list($key,$value) = $pair;
				if(!in_array($key,$all_keys)) $all_keys[] = $key;
			}
		}
		++$row_count;
	}
	
	$all_keys[] = "flags"; // Use key flags for flags that can occur in INFO column
	
	
	//2D array that contains key=>value pairs, "." if empty
	$all_cols = array_fill_keys($all_keys,array_fill(0,$row_count,"."));
	
	for($i = 0; $i<$row_count; ++$i)
	{
		$parts = explode(";", $info_col[$i]);
		$flags = array();
		foreach($parts as $raw_part)
		{
			$pair = explode("=", $raw_part);
			
			$key = "";
			$value = "";
			
			if(count($pair) == 2)
			{
				list($key,$value) = $pair;
			}
			else
			{
				$flags[] = $pair[0];
			}
			
			if(array_key_exists($key,$all_cols))
			{
				$all_cols[$key][$i] = $value;
			}
			
		}
		if(count($flags) > 0) $all_cols["FLAGS"][$i] = implode(";",$flags);
	}
	
	return $all_cols;
}

################################## MAIN ##################################
$bedpe_in = Matrix::fromTSV($in);

$headers = $bedpe_in->getHeaders();


//parse headers
//Remove comments that describe "contigs" since they only describe length of chromosomes from reference genome
$bedpe_in->removeComment("#contig=", true);
//Replace fileformat header line
$bedpe_in->removeComment("#fileformat=BEDPE");
$bedpe_in->prependComment("#fileformat=BEDPE_TUMOR_NORMAL_PAIR");
$bedpe_in->prependComment("#ANALYSISTYPE=MANTA_TUMOR_NORMAL_PAIR");


################################## Parse columns from INFO_A column ##################################
//Expand certain information for INFO_A colum
//We skip INFO_B because it includes only infos in case of BNDs, and these infos are almost the same as in INFO_A in that rare case
$info_fields = get_info_fields($bedpe_in->getComments());
$i_info_a = $bedpe_in->getColumnIndex("INFO_A");

$results = array("SOMATICSCORE" => array(), "JUNCTION_SOMATICSCORE" => array(), "FLAGS" => array());
if($bedpe_in->rows() > 0) $results = (expand_info_col($bedpe_in->getCol($i_info_a)));

//Add information for SOMATICSCORE, JUNCTION_SOMATICSCORE and FLAGS to output file.
if(array_key_exists("SOMATICSCORE",$results))
{
	$bedpe_in->addCol($results["SOMATICSCORE"], "SOMATICSCORE", $info_fields["SOMATICSCORE"]["Description"]);
}
else //if nothing set in file, add empty SOMATICSCORE column
{
	$bedpe_in->addCol(array_fill(0,$bedpe_in->rows(), "."), "SOMATICSCORE", "Somatic variant quality score");
}
if(array_key_exists("JUNCTION_SOMATICSCORE",$results))
{
	$bedpe_in->addCol($results["JUNCTION_SOMATICSCORE"], "JUNCTION_SOMATICSCORE", $info_fields["JUNCTION_SOMATICSCORE"]["Description"]);
}
else //if nothing set in file, add empty JUNCTION_SOMATICSCORE column
{
	$bedpe_in->addCol(array_fill(0,$bedpe_in->rows(), "."), "JUNCTION_SOMATICSCORE", "If the SV junction is part of an EVENT (ie. a multi-adjacency variant), this field provides the SOMATICSCORE value for the adjacency in question only.");
}
if(array_key_exists("FLAGS", $results))
{
	$bedpe_in->addCol($results["FLAGS"], "FLAGS", "Flags that occur in INFO column.");
}


##################################PARSE columns from FORMAT column ##################################
$i_format = $bedpe_in->getColumnIndex("FORMAT",false,true);
$i_tumor = $bedpe_in->getColumnIndex($tid, false, true);

$cols = array("PR_REF" => array(), "SR_REF" => array(), "PR_ALT" => array(), "SR_ALT" => array());
if($bedpe_in->rows() > 0)
{
	$cols = expand_format_col( $bedpe_in->getCol($i_format), $bedpe_in->getCol($i_tumor) );
}

$bedpe_in->addCol($cols["PR_REF"], "TUM_PR_REF", "Spanning paired-read support for the ref alleles in {$tid}.");
$bedpe_in->addCol($cols["SR_REF"], "TUM_SR_REF", "Split reads for the ref alleles in the order listed, for reads where P(allele|read)>0.999 in {$tid}.");
$bedpe_in->addCol($cols["PR_ALT"], "TUM_PR_ALT", "Spanning paired-read support for the alt alleles in {$tid}.");
$bedpe_in->addCol($cols["SR_ALT"], "TUM_SR_ALT", "Split reads for the alt alleles in the order listed, for reads where P(allele|read)>0.999 in {$tid}.");

//Annotate supporting Alleles of tumor sample
$i_normal = $bedpe_in->getColumnIndex($nid, false, true);
if($bedpe_in->rows() > 0)
{
	$cols = expand_format_col( $bedpe_in->getCol($i_format), $bedpe_in->getCol($i_normal) );
}

$bedpe_in->addCol($cols["PR_REF"], "NOR_PR_REF", "Spanning paired-read support for the ref alleles in {$nid}.");
$bedpe_in->addCol($cols["PR_ALT"], "NOR_PR_ALT", "Spanning paired-read support for the alt alleles in {$nid}.");
$bedpe_in->addCol($cols["SR_REF"], "NOR_SR_REF", "Split reads for the ref alleles in the order listed, for reads where P(allele|read)>0.999 in {$nid}.");
$bedpe_in->addCol($cols["SR_ALT"], "NOR_SR_ALT", "Split reads for the alt alleles in the order listed, for reads where P(allele|read)>0.999 in {$nid}.");

//Remove original FORMAT and corresponding data columns, including comments
$bedpe_in->removeCol($i_tumor);
$bedpe_in->removeCol($i_normal);
$bedpe_in->removeCol($i_format);
$bedpe_in->removeComment("#FORMAT=<ID=SR", true);
$bedpe_in->removeComment("#FORMAT=<ID=PR", true);

//Remove cols for + and - strand because they do not make much sense in NGS sequencing
$bedpe_in->removeCol($bedpe_in->getColumnIndex("STRAND_A"));
$bedpe_in->removeCol($bedpe_in->getColumnIndex("STRAND_B"));


################################## Reorder Columns ##################################
//Move unimportant columns to the end
$to_be_moved = array("FLAGS","REF_A","ALT_A","REF_B","ALT_B","INFO_A","INFO_B","NAME_A","NAME_B");
foreach($to_be_moved as $idx)
{
	$bedpe_in->moveColToEnd($bedpe_in->getColumnIndex($idx));
}

//Remove ID column because ID is included in NAME_A or NAME_B column already
$bedpe_in->removeCol($bedpe_in->getColumnIndex("ID"));
//Remove "QUAL" column if empty (we have SOMATICSCORE and JUNCTION_SOMATICSCORE for tumor-normal-pairs
$qual_is_empty = true;
$i_qual = $bedpe_in->getColumnIndex("QUAL");
for($i=0; $i<$bedpe_in->rows(); ++$i)
{
	if($bedpe_in->get($i,$i_qual) != "." && !empty($bedpe_in->get($i,$i_qual)))
	{
		$qual_is_empty = false;
		break;
	}
}
if($qual_is_empty)
{
	$bedpe_in->removeCol($i_qual);
}


$bedpe_in->toTSV($out);

?>