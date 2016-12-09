<?php

/** 
	@page merlin_generate_model
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$parser = new ToolBase("merlin_generate_model", "\$Rev: 2 $", "Generates model-files to use with merlin's parametric analysis.");
$parser->addInfile("ped",  "Input PED file.", false);
$parser->addInfile("dat",  "Input DAT file", false);
$parser->addEnum("type", "Linkage model type.", false, array("dominant_hp"));
$parser->addOutfile("out",  "Output MODEL file", true);
extract($parser->parse($argv));


if ($type=="dominant_hp")
{
	//calculate desease allel frequency, assuming affected have Aa, not affected aa
	$pedfile = file($ped);
	$datfile = file($dat);
	
	$affection_marker_pos=5;	//the first 5 entries in a ped file are no markers
	$affection_found=False;
	foreach($datfile as $line)
	{
		$exploded_line=explode("\t",$line);
		$affection_found = (trim($exploded_line[1])=="affection_status");
		if ($affection_found) break;
		
		$affection_marker_pos++;
	}		
	
	if (!$affection_found)
	{
		trigger_error("No affection status found in $dat", E_USER_ERROR);
	}
	
	$disease_allels=0;
	$wildtype_allels=0;
	
	foreach($pedfile as $line)
	{		
		$exploded_line=explode("\t",$line);
		if (count($exploded_line)<$affection_marker_pos) continue;
		if ($exploded_line[$affection_marker_pos]==2) 		//affected, assume genotype Aa
		{
			$disease_allels++;
			$wildtype_allels++;
		}
		else if ($exploded_line[$affection_marker_pos]==1)	//not affected, assume genotype aa
		{
			$wildtype_allels++;
			$wildtype_allels++;
		}
	}
	
	$allel_freq=$disease_allels/($wildtype_allels+$disease_allels);
	file_put_contents($out, implode ("\t", array("affection_status",   $allel_freq, "0.0001,1.0,1.0", "Dominant_high_pen")));
}
else 
{
	trigger_error("model type not supported", E_USER_ERROR);
}

?>