<?php

/**
	@page chip_ngs_compare
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("chip_ngs_compare", "\$Rev: 712 $", "Array and NGS SNP comparison.");
$parser->addInfile("chip",  "Input SNP array genotype file in TSV format.", false);
$parser->addInfile("ngs",  "Input NGS variant list in TSV format.", false);
$parser->addEnum("type", "Chip type.", false, array("illumina_6k", "affymetrix_6.0", "affymetrix_250k", "cytoscan_hd"));
//optional
$parser->addOutfile("out",  "Output text file (stdout if unset).", true);
extract($parser->parse($argv));

//load affymetrix meta data - mapping of affy SNP names to RS numbers, chromosomal location
/*
#name   rs      chr     pos     cM
SNP_A-1780419   rs6576700       1       84875173        102.95735950781
SNP_A-1780418   rs17054099      5       156390980       154.762533280549
SNP_A-1780415   rs7730126       5       158729947       156.668365937775
*/
function load_affy($filename)
{
	$output = array();
	
	$file = file($filename);
	foreach($file as $line)
	{
		//skip non-content lines
		if (trim($line)=="" || $line[0]=="#") continue;
		
		//skip SNPs without position information
		list($name, , $chr, $pos, , $a, $b) = explode("\t", $line);
		if ($chr=="---") continue;
		
		$output[$name] = array(trim($chr), trim($pos), trim($a), trim($b));
	}
	return $output;
}

function calculate_counts($title, $genotypes, &$output)
{
	$count = count($genotypes);
	
	//count het
	$geno_counts = array_count_values($genotypes);
	$het = 0;
	foreach($geno_counts as $geno => $geno_count)
	{
		if ($geno[0]!=$geno[1]) $het += $geno_count;
	}
	
	//general info
	$output[] = "$title SNPs: $count (".number_format(100*$het/$count,2)."% het)\n";
	
	if (!preg_match("/^[ACGT]{2}$/", $geno))
	{
		trigger_error("Invalid genotype '$geno' found!", E_USER_ERROR);
	}
}

$chip_data = array();
if ($type=="affymetrix_6.0")
{
	$chip_data = load_affy(get_path("data_folder")."/dbs/affymetrix/GenomeWideSNP_6.na32.annot.tsv.new");
}
else if ($type=="affymetrix_250k")
{
	$chip_data = load_affy(get_path("data_folder")."/dbs/affymetrix/Mapping250K_Nsp.na32.annot.tsv.new");
}
else if ($type=="illumina_6k")
{
	$chip_data = load_affy(get_path("data_folder")."/dbs/illumina/6k_meta.tsv.new");
}
else if ($type=="cytoscan_hd")
{
	$chip_data = load_affy(get_path("data_folder")."/dbs/affymetrix/CytoScanHD_Array.na32.1.annot.tsv.new");
}

$output = array();

//load chip genotypes
/*
Probe Set ID	112445_(GenomeWideSNP_6).birdseed-v2.chp Call Codes	112445_(GenomeWideSNP_6).birdseed-v2.chp Confidence
SNP_A-2131660	AB	0.018
SNP_A-1967418	BB	0.009
*/
$chp_geno = array();
$file = file($chip);
foreach($file as $line)
{
	//skip non-content lines
	if (trim($line)=="" || $line[0]=="#") continue;
	
	//skip SNPs without genotype call
	list($id, $geno) = explode("\t", $line);
	if ($geno=="NoCall") continue;
	
	//skip SNPs with no position info
	if (!isset($chip_data[$id])) continue;
	
	//determine genotype
	list($chr, $pos, $a, $b) = $chip_data[$id];

	$chp_geno["chr".$chr."_".$pos] = strtr($geno, array("AA"=>"$a$a", "AB"=>"$a$b", "BB"=>"$b$b"));
}
calculate_counts("CHP", $chp_geno, $output);

//load NGS genotypes
$ngs_geno = array();
$file = file($ngs);
foreach($file as $line)
{
	//skip non-content lines
	if (trim($line)=="" || $line[0]=="#") continue;
	
	//determine genotype
	list($chr, $pos, , $ref, $obs, $geno, , , , , , , , , , , $rs) = explode("\t", $line);
	
	if ($ref=="-" || $obs=="-" || $rs==".") continue;
	
	$ngs_geno[$chr."_".$pos] = strtr($geno, array("hom"=>"$obs$obs", "het"=>"$ref$obs"));
}
calculate_counts("NGS", $ngs_geno, $output);

//count matches
$matches = 0;
$same = 0;
foreach($ngs_geno as $pos => $geno)
{
	if (isset($chp_geno[$pos]))
	{
		++$matches;
		if ($geno==$chp_geno[$pos])
		{
			++$same;
		}
	}
}
$output[] = "OVERLAP: $matches\n";
$output[] = "SAME GENOTYPE: $same (".number_format(100*$same/$matches,2)."%)\n";

if (!isset($out))
{
	print implode("", $output);
}
else
{
	file_put_contents($out, $output);
}

?>
