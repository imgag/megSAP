<?php

/**
	@page chip2plink
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("chip2plink", "\$Rev: 712 $", "Converts a SNP array genotype file to PLINK input files (MAP and PED) using meta data from the mapping array table.");
$parser->addInfile("chip",  "SNP array genotype file.", false);
$parser->addInfile("meta",  "Text export of mapping array table.", false);
$parser->addEnum("type", "Chip type.", false, array("illumina_6k", "affymetrix_6.0", "affymetrix_250k", "cytoscan_hd"));
$parser->addString("family", "Family identifier string.", false);
$parser->addOutfile("out_map", "Output mapping file name.", false);
$parser->addOutfile("out_ped", "Output pedigree file name.", false);
$parser->addOutfile("out_fam", "Text file with estimated family relations.", false);
$parser->addFlag("auto_only", "Writes only autosomal SNPs to output files.");
extract($parser->parse($argv));

//load affymetrix meta data - mapping of affy SNP names to RS numbers, chromosomal location and cM value (if necessary)
function load_affy($filename)
{
	$output = array();
	
	$file = file($filename);
	foreach($file as $line)
	{
		list($name, $rs, $chr, $pos, $cm) = explode("\t", $line);
		$output[$name] = array(trim($rs), trim($chr), trim($pos), trim($cm));
	}
	return $output;
}

$chip_data = array();
if ($type=="affymetrix_6.0")
{
	$chip_data = load_affy(get_path("data_folder")."/dbs/affymetrix/GenomeWideSNP_6.na32.annot.tsv");
}
else if ($type=="affymetrix_250k")
{
	$chip_data = load_affy(get_path("data_folder")."/dbs/affymetrix/Mapping250K_Nsp.na32.annot.tsv");
}
else if ($type=="illumina_6k")
{
	$chip_data = load_affy(get_path("data_folder")."/dbs/illumina/6k_meta.tsv");
}
else if ($type=="cytoscan_hd")
{
	$chip_data = load_affy(get_path("data_folder")."/dbs/affymetrix/CytoScanHD_Array.na32.1.annot.tsv");
}

// write map file
$chp = Matrix::fromTSV($chip);
$map = new Matrix();
$skipped = array();
$num = $chp->rows();
for ($i=1; $i<$num; ++$i)
{
	$row = $chp->getRow($i);
	if ($type=="illumina_6k")
	{
		list(, $name, , $chr, $pos) = $row;
		$CM = "---";
		if (isset($chip_data[$name]))
		{
			$cM = $chip_data[$name][3];
		}
	}
	else if ($type=="affymetrix_6.0" || $type=="affymetrix_250k" || $type=="cytoscan_hd")
	{
		if (!isset($chip_data[$row[0]]))
		{
			$skipped[$row[0]] = true;
			continue;
		}
		
		list($name, $chr, $pos, $cM) = $chip_data[$row[0]];
	}
	
	if (!$auto_only || ctype_digit($chr))
	{
		$map->addRow(array($chr, $name, $cM, $pos));
	}
}
if (count($skipped)!=0)
{
	trigger_error(count($skipped)." of $num SNPs not found in Affymetrix annotation file!", E_USER_WARNING);
}
$map->toTSV($out_map);


//create pedigree file (base info) and look up person indices in chip file
$data = Matrix::fromTSV($meta);
$ped = new Matrix();
$headers = $chp->getRow(0);
$indices = array();
$num = $data->rows();
for ($i=0; $i<$num; ++$i)
{
	$row = $data->getRow($i);
	$fam_id = $row[6];
	
	if ($fam_id == $family)
	{
		$per_id = $row[1];
		$father = trim($row[10]);
		if ($father=="") $father="0";
		$mother = trim($row[9]);
		if ($mother=="") $mother="0";
		$valid_genders = array("m"=>"1", "w"=>"2", "unbekannt"=>"other");
		$gender = strtr(trim($row[7]), $valid_genders);
		if (!in_array($gender, array_values($valid_genders)))
		{
			trigger_error("Unknown gender '$gender'!", E_USER_ERROR);
		}
		$valid_status = array("erkrankt"=>"2", "gesund"=>"1", "unbekannt"=>"0");
		$affected = strtr(trim($row[8]), $valid_status);
		if (!in_array($affected, array_values($valid_status)))
		{
			trigger_error("Unknown status '$affected'!", E_USER_ERROR);
		}
		$ped->addRow(array(strtr($fam_id, array(" "=>"_")), $per_id, $father, $mother, $gender, $affected));
		
		// look up patient index in chip file
		$mft_id = $row[12];
		$index = array();
		for($j=0; $j<count($headers); ++$j)
		{
			if ($type=="illumina_6k" && contains($headers[$j], $mft_id.".GType"))
			{
				$index[] = $j;
			}
			else if ($type=="affymetrix_6.0"  && starts_with($headers[$j], $mft_id."_") && contains($headers[$j], "Call") && !contains($headers[$j], "Forward Strand"))
			{
				$index[] = $j;
			}
			else if ($type=="affymetrix_250k" && starts_with($headers[$j], $mft_id."_") && contains($headers[$j], "Call") && !contains($headers[$j], "Forward Strand"))
			{
				$index[] = $j;
			}
			else if ($type=="cytoscan_hd" && (starts_with($headers[$j], "Call Codes") || (starts_with($headers[$j], $mft_id."_") && contains($headers[$j], "Call Codes"))))
			{
				$index[] = $j;
			}
		}
		$count = count($index);
		if ($count!=1)
		{
			trigger_error("MFT id '$mft_id' found '$count' times in chip file - expected once!", E_USER_ERROR);
		}
		$indices[] = $index[0];
	}
}

// abort if noone of that family is found
if (count($indices)==0)
{
	trigger_error("Family '$family' not found in meta data!", E_USER_ERROR);
}

// extend pedigree file with genotypes and write it
$replacements = array("AA"=>"1\t1", "BB"=>"2\t2", "AB"=>"1\t2", "BA"=>"2\t1", "NC"=>"0\t0", "NoCall"=>"0\t0");
$num = $chp->rows();
for ($i=1; $i<$num; ++$i)
{
	$row = $chp->getRow($i);
	
	//skip SNPs with no annotation
	if (isset($skipped[$row[0]]))
	{
		continue;
	}
	
	//skip non-autosomal SNPs if auto_only is activated
	if ($auto_only)
	{
		if ($type=="illumina_6k" && !ctype_digit($row[3]))
		{
			continue;
		}
		else if (($type=="affymetrix_250k" || $type=="affymetrix_6.0" || $type=="cytoscan_hd") && !ctype_digit($chip_data[$row[0]][1]))
		{
			continue;
		}
	}

	$ped_col = array();
	foreach($indices as $index)
	{
		$ped_col[] = strtr($row[$index], $replacements);
	}
	$ped->addCol($ped_col);
}
$ped->toTSV($out_ped);

// returns the fraction of heterocygous genotypes
function het_fraction($col, $auto_indices)
{
	$num = count($auto_indices);
	$count_het = 0;
	for ($i=0; $i<$num; ++$i)
	{
		$index = $auto_indices[$i];
		if ($col[$index]=="1\t2" || $col[$index]=="2\t1")
		{
			++$count_het;
		}
	}
	
	return $count_het / $num;
}

// check family information based on correlation of genotypes and medelian error fraction
function relation($col1, $col2, $auto_indices)
{
	$replacements = array("1\t1"=>1, "2\t2"=>3, "1\t2"=>2, "2\t1"=>2);

	$c1 = array();
	$c2 = array();
	$num = count($auto_indices);
	for ($i=0; $i<$num; ++$i)
	{
		$index = $auto_indices[$i];
		if ($col1[$index]=="0\t0" || $col2[$index]=="0\t0")
		{
			continue;
		}
		
		$c1[] = strtr($col1[$index], $replacements);
		$c2[] = strtr($col2[$index], $replacements);		
	}
	
	//correlation
	$corr = correlation($c1, $c2);
	
	//mendelian error fraction
	$me = 0; 
	for($k=0; $k<count($c1); ++$k)
	{
		if (($c1[$k]==1 && $c2[$k]==3) || ($c1[$k]==3 && $c2[$k]==1))
		{
			++$me;
		}
	}
	$mef = $me/count($c1);
	
	//not related (2.5 sigma above/below mean)
	if ($corr<0.27 || $mef>0.089)
	{
		return array("-", $corr, $mef);
	}
	
	//direct relation (parent-child)
	if ($corr>0.5 && $mef<0.01)
	{
		return array("parent-child", $corr, $mef);
	}
	
	//siblings
	if ($corr>0.5 && $mef<0.04)
	{
		return array("siblings", $corr, $mef);
	}
	
	//distant relation
	return array("distant", $corr, $mef);
}

// get array of x chromosomal genotypes
function x_genotypes($col, $x_indices)
{
	$replacements = array("1\t1"=>"AA", "2\t2"=>"BB", "1\t2"=>"AB", "2\t1"=>"AB");

	$output = array();
	for ($i=0; $i<count($x_indices); ++$i)
	{
		$index = $x_indices[$i];
		if ($col[$index]=="0\t0" || $col[$index]=="0\t0")
		{
			continue;
		}
		
		$output[] = strtr($col[$index], $replacements);
	}
	
	return $output;
}


//look up:
// - indices of autosomal SNPs - for family relation check
// - indices of X SNPs - for gender check
$auto_indices = array();
$x_indices = array();
for ($i=0; $i<$map->rows(); ++$i)
{
	$row = $map->getRow($i);
	if (ctype_digit($row[0]))
	{
		$auto_indices[] = $i;
	}
	if ($row[0]=="X")
	{
		$x_indices[] = $i;
	}
}

//verify genders given in meta data file
$details = array();
$details[] = "#IID	gender	het_ratio_x	het_ratio_auto	affected";
$messages = array();
$num = $ped->rows();
for ($i=0; $i<$num; ++$i)
{
	$row = $ped->getRow($i);
	$info = array_splice($row, 0, 6);
	$iid = $info[1];
	$ref = $info[4];
	$affected = $info[5]-1;
	list($obs, $ratio_x) = gender(x_genotypes($row, $x_indices), "AB", 0.12, 0.12);
	$details[] = "$iid	$obs	".number_format($ratio_x,4)."	".number_format(het_fraction($row, $auto_indices),4)."	$affected";
	
	if ($ref=="other")
	{
		$messages[] = "#Notice: Gender of '$iid' not given. It is calculated as '$obs'!";
	}
	else if ($ref==2 && $obs=='m')
	{
		$messages[] = "#Warning: Gender of '$iid' is given as '$ref' but calculated as '$obs'!";
	}
	else if ($ref==1 && $obs=='f')
	{	
		$messages[] = "#Warning: Gender of '$iid' is given as '$ref' but calculated as '$obs'!";
	}									
}

//verify relations given in meta data file
$fam = array();
$details[] = "#IID1	IID2	relation	correlation	mendelian_error_fraction";
for ($i=0; $i<$num; ++$i)
{
	$row1 = $ped->getRow($i);
	$info1 = array_splice($row1, 0, 6);
	
	for ($j=$i+1; $j<$num; ++$j)
	{
		$row2 = $ped->getRow($j);
		$info2 = array_splice($row2, 0, 6);
		
		list($relation, $corr, $mef) = relation($row1, $row2, $auto_indices);
		$details[] = $info1[1]."	".$info2[1]."	$relation	".number_format($corr,4)."	".number_format($mef,4);
		
		$fam[$info1[1]][$info2[1]] = $relation;
		
		// unrelated (parents)
		if ($info1[2] == "0" && $info1[3] == "0" && $info2[2] == "0" && $info2[3] == "0")
		{
			if ($relation!="-")
			{
				$messages[] = "#Warning: Unrelated persons ".$info1[1]." and ".$info2[1]." classified as '$relation'!";
			}
		}
		// siblings
		else if ($info1[2] == $info2[2] && $info1[3] == $info2[3])
		{
			if ($relation!="siblings")
			{
				$messages[] = "#Warning: Siblings ".$info1[1]." and ".$info2[1]." classified as '$relation'!";
			}		
		}
		// parent-child
		else if ($info1[1] == $info2[2] || $info1[1] == $info2[3] || $info2[1] == $info1[2] || $info2[1] == $info1[3])
		{
			if ($relation!="parent-child")
			{
				$messages[] = "#Warning: Parent-child ".$info1[1]." and ".$info2[1]." classified as '$relation'!";
			}
		}
		else
		{
			$messages[] = "#Warning: Unknown family relation between ".$info1[1]." and ".$info2[1]."!";
		}
	}
}
$details[] = "";

file_put_contents($out_fam, implode("\n", array_merge($messages, $details)));

?>
