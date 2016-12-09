<?php

/**
	@page linkage

*/

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";

require_once($basedir."Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("linkage", "\$Rev: 2 $", "Linkage pipeline.");
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
	$parser->execTool("php ".$basedir."Chips/filter_confidence.php", "-in $geno -thres $thres -out $geno_filtered");
}
else
{
	$parser->exec("cp", "$geno $geno_filtered", false);
}

// create map/ped file
$fam = $family."_family.txt";
$parser->execTool("php ".$basedir."Chips/chip2plink.php", "-chip $geno_filtered -meta $meta -type $type -family $family -out_map merlin.map -out_ped merlin.ped -out_fam $fam");

// remove last column from map file
$file = Matrix::fromTSV("merlin.map");
$file->removeCol(3);
$file->toTSV("merlin.map");

//create dat file
$output = array();
$output[] = "A	affection_status";
$file = Matrix::fromTSV("merlin.map");
for ($i=0; $i<$file->rows(); ++$i)
{
	$row = $file->getRow($i);
	$output[] = "M	".$row[1];
}
file_put_contents("merlin.dat", implode("\n", $output)); 

//find unlikely genotypes using merlin, save them in merlin.err (filename hard coded in merlin, no other filename possible),
$parser->exec(get_path("merlin"), "merlin.dat -p merlin.ped -m merlin.map --error", true);

//remove unkikely genotypes from analysis using pedwipe, save results in wiped.dat and wiped.ped (filenames hard coded in pedwipe, no other filename possible)
$parser->exec(get_path("pedwipe"), "-d merlin.dat -p merlin.ped", false);

//generate model file for linkage analysis
$parser->execTool("php ".$basedir."Chips/merlin_generate_model.php", "-type dominant_hp -dat merlin.dat -ped merlin.ped -out merlin.model");

// apply merlin to calculate linkage in a paramertic model
$parser->exec(get_path("merlin"), "-d wiped.dat -p wiped.ped -m merlin.map --model merlin.model  --markerNames --pdf --tabulate --prefix ".$family." --perFamily", true);

//build bedfile based on analysis results
$parser->execTool("php ".$basedir."Chips/merlin_report.php", "-in ".$family."-parametric.tbl -out _region.bed");

//build annotated file based on bedfile
$parser->execTool("php ".$basedir."Chips/bed_annotation.php", "-in _region.bed -out merlin.tsv");

//create and store report
$output = array();
$output[] = "Report zur Linkageanalyse:";
$output[] = "";
$output[] = "Familie: $family";
$output[] = "Datum: ".date("d.m.Y");
$output[] = "Revision der Analysepipeline: ".get_svn_rev();
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

file_put_contents($family."_report.txt", implode("\n", $output)); 

?>
