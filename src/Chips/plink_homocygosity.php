<?php

/**
	@page plink_homocygosity
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("plink_homocygosity", "\$Rev: 877 $", "Wrapper for PLINK homocygosity mapping and converter to improved output format.");
$parser->addInfile("in_map",  "Input mapping file name.", false);
$parser->addInfile("in_ped",  "Input pedigree file name.", false);
$parser->addOutfile("out", "Output homocygosity mapping result file name.", false);
$parser->addInt("w_snp", "Size of the window (in SNPs).", false);
$parser->addInt("w_kb", "Size of the window (in kb).", false);
$parser->addInt("w_het", "Number of heterozygotes allowed in a window.", false);
$parser->addInt("w_miss", "Number of missing genotype calls allowed in the window.", false);
$parser->addFloat("w_thres", "Ratio of windows that must be homocygote to consider a base homocygote.", false);
$parser->addInt("s_snp", "Mimum size of segments (in SNPs).", false);
$parser->addInt("s_kb", "Mimum size of segments (in kb).", false);
$parser->addInt("s_dens", "Required minimum density of segments (in 1 SNP per x kb).", false);
$parser->addInt("s_gap", "Maximum size of gaps in segments (in kb).", false);
extract($parser->parse($argv));

// copy to temporary file location
$temp = $parser->tempFile();
copy2($in_map, $temp.".map");
copy2($in_ped, $temp.".ped");

// execute plink
$plink = get_path("plink");
$parser->exec($plink, "--file $temp --homozyg --noweb --homozyg-window-kb $w_kb --homozyg-window-snp $w_snp --homozyg-window-het $w_het --homozyg-window-missing $w_miss --homozyg-window-threshold $w_thres --homozyg-snp $s_snp --homozyg-kb $s_kb --homozyg-density $s_dens --homozyg-gap $s_gap --out $temp", false);

// load PLINK output file and fix format errors:
// - column format with multiples of four spaces
// - chr23 should be chrX
// - phenotypes 1.000/2.000 > 1/2
$output = new Matrix();
$file = file($temp.".hom");
for ($i=0; $i<count($file); ++$i)
{
	$line = $file[$i];
	$line = preg_replace("/\s\s+/", " ", " ".$line);
	
	$parts = explode(" ", nl_trim($line));
	array_splice($parts, 0, 1);
	
	$parts[3] = strtr($parts[3], array("23"=>"X"));
	$parts[2] = strtr($parts[2], array("1.000"=>"1", "2.000"=>"2"));
	
	if ($i==0)
	{
		$output->setHeaders($parts);
	}
	else
	{
		$output->addRow($parts);	
	}
}

// store fixed format
$output->toTSV($out);

?>
