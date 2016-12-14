<?php

/**
	@page homocygosity

*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("homocygosity", "Homocygosity mapping pipeline.");
$parser->addInfile("geno",  "Genoptypes input file.", false);
$parser->addInfile("meta",  "Meta data file.", false);
$parser->addEnum("type", "Chip type.", false, array("illumina_6k", "affymetrix_6.0", "affymetrix_250k", "cytoscan_hd"));
$parser->addString("family", "Family identifier string.", false);
$parser->addFloat("thres", "Confidence value threshold (only used for Affymetrix 6.0 array).", true, 0.02);
extract($parser->parse($argv));


// confidence value filter
$geno_filtered = $family."_geno_filtered.txt";
if ($type=="affymetrix_6.0")
{
	$parser->execTool("Chips/filter_confidence.php", "-in $geno -thres $thres -out $geno_filtered");
}
else
{
	$parser->exec("cp", "$geno $geno_filtered", false);
}

// convert to PLINK format
$map = $family."_map.txt";
$ped = $family."_ped.txt";
$fam = $family."_family.txt";
$parser->execTool("Chips/chip2plink.php", "-chip $geno_filtered -meta $meta -type $type -family $family -out_map $map -out_ped $ped -out_fam $fam");

// apply plink
$params = "";
if ($type=="affymetrix_6.0" || $type=="cytoscan_hd")
{
	$params = "-w_snp 60 -w_kb 2000 -w_het 8 -w_miss 10 -w_thres 0.5 -s_snp 10 -s_kb 1000 -s_dens 100  -s_gap 5000";
}
else if ($type=="affymetrix_250k")
{
	$params = "-w_snp 40 -w_kb 1500 -w_het 4 -w_miss 7 -w_thres 0.5 -s_snp 10 -s_kb 1000 -s_dens 500 -s_gap 10000";
}
else if ($type=="illumina_6k")
{
	$params = "-w_snp 10 -w_kb 1000 -w_het 1 -w_miss 2 -w_thres 0.5 -s_snp 10 -s_kb 1000 -s_dens 1000 -s_gap 20000";
}
else
{
	trigger_error("No PLINK parameters known for chip type '$type'. Aborting!", E_USER_ERROR);
}


$parser->addInt("w_snp", "Size of the window (in SNPs).", false);
$parser->addInt("w_kb", "Size of the window (in kb).", false);
$parser->addInt("w_het", "Number of heterozygotes allowed in a window.", false);
$parser->addInt("w_miss", "Number of missing genotype calls allowed in the window.", false);
$parser->addFloat("w_thres", "Ratio of windows that must be homocygote to consider a base homocygote.", false);
$parser->addInt("s_snp", "Mimum size of segments (in SNPs).", false);
$parser->addInt("s_kb", "Mimum size of segments (in kb).", false);
$parser->addInt("s_dens", "Required minimum density of segments (in 1 SNP per x kb).", false);
$parser->addInt("s_gap", "Maximum size of gaps in segments (in kb).", false);

$hom = $family."_hom.txt";
$parser->execTool("Chips/plink_homocygosity.php", "-in_map $map -in_ped $ped -out $hom $params");


//create diagrams
$parser->execTool("Chips/plink_diagrams.php", "-in $hom");

//create regions
$reg = $family."_regions.bed";
$parser->execTool("Chips/plink_regions.php", "-in $hom -out $reg");

//annotate regions
$anno = $family."_annotation.tsv";
$parser->execTool("Chips/bed_annotation.php", "-in $reg -out $anno");

//create and store report
$output = array();
$output[] = "Report zur Homozygotieanalyse:";
$output[] = "";
$output[] = "Familie: $family";
$output[] = "Datum: ".date("d.m.Y");
$output[] = "Revision der Analysepipeline: ".repository_revision();
$output[] = "Chip: $type";
$output[] = "";
$output[] = "Patienten:";
$output[] = "  #Name (Geschlecht, Status)	ID-Patient	ID-Vater	ID-Mutter";
$meta_file = file($meta);
foreach($meta_file as $line)
{
	list(, $dna, , $name, $prename, , $family_string, $gender, $status, $mother, $father) = explode("\t", trim($line));
	if ($family_string != $family) continue;
	if ($mother=="") $mother = "-";
	if ($father=="") $father = "-";
	$output[] = "  $prename $name ($gender, $status)	$dna	$father	$mother";
}
$output[] = "";
$output[] = "Homozygote Intervalle (hg19):";
$output[] = "  #Chromosom	Start	Ende	MB";
$anno_file = file($anno);
$inter_done = array();
foreach($anno_file as $line)
{
	if ($line[0]=='#') continue;
	list($inter) = explode("\t", trim($line));
	if (!in_array($inter, $inter_done))
	{
		list($chr, $start, $end) = explode("\t", strtr($inter, "-:", "\t\t"));
		$output[] = "  $chr	$start	$end	".number_format(($end-$start)/1000000, 2);
		$inter_done[] = $inter;
	}
}
$output[] = "";

file_put_contents($family."_report.txt", implode("\n", $output)); 

?>
