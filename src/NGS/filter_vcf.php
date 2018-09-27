<?php
/**
	@page filter_vcf
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("filter_vcf", "Flag variants in a VCF file according to quality.");
$parser->addInfile("in", "Input variant file in VCF format.", false);
$parser->addString("t_id",  "Tumor ID in VCF file.", false);
$parser->addString("n_id",  "Normal ID in VCF file.", false);
$parser->addOutfile("out", "Output variant file in VCF format.", false);
extract($parser->parse($argv));

$in_file = Matrix::fromTSV($in);

//check filter column present 
$fc = $in_file->getColumnIndex("FILTER", false, false);
if($fc!=6)	trigger_error("Wrong column index for filter column! Is ".(empty($fc)?"empty":$fc)." should be 6.", E_USER_ERROR);

//determine tumor/normal column
$tumor_col = $in_file->getColumnIndex($t_id);
if(is_null($tumor_col)) trigger_error("Could not identify tumor column in $in.", E_USER_ERROR);
$normal_col = $in_file->getColumnIndex($n_id);
if(is_null($normal_col)) trigger_error("Could not identify normal column in $in.", E_USER_ERROR);

//process variants
for($i=0;$i<$in_file->rows();++$i)
{
	//init
	$row = $in_file->getRow($i);
	$chr = $row[0];
	$start = $row[1];
	$end = $start;
	$ref = $row[3];
	$alt = $row[4];
	$genotype = $row[8];
	
	//determine variant type
	$type = "SNV";
	if(strlen($ref) > 1)	$type = "INDEL";
	foreach(explode(",", $alt) as $a)
	{
		if(strlen($a)>1) $type = "INDEL";
	}
	
	//correct indels
	if($type == "INDEL" && $ref!="-" && $alt!="-")	list($start, $end, $ref, $alt) = correct_indel($start, $ref, $alt);
	
	//extract depth/frequency
	$td = NULL;
	$tf = NULL;
	$nd = NULL;
	$nf = NULL;
	if($type == "SNV" && preg_match("/^[acgtACGT]*$/", $alt))
	{
		list($td, $tf) = vcf_strelka_snv($genotype, $row[$tumor_col], $alt);
		list($nd, $nf) = vcf_strelka_snv($genotype, $row[$normal_col], $alt);
	}
	else if($type == "INDEL")
	{
		list($td, $tf) = vcf_strelka_indel($genotype, $row[$tumor_col]);
		list($nd, $nf) = vcf_strelka_indel($genotype, $row[$normal_col]);
	}
	//print "$chr:$start $ref>$alt $type T=$td/$tf N=$nd/$nf\n";
	
	//filter by depth/frequency (if possible)
	$filters = array();
	if ($td!=NULL && $tf!=NULL && $nd!=NULL && $nf!=NULL)
	{
		if ($td*$tf<2.9) $filters[] = "lt-3-reads";
		if ($td<20) $filters[] = "depth-tum";
		if ($nd<20) $filters[] = "depth-nor";
		if($tf<0.05) $filters[] = "freq-tum";
		if($nf>1/6*$tf) $filters[] = "freq-nor";
	}
	
	//extract old filter entries
	$filters_old = explode(";", trim($row[6]));
	if(count($filters_old)==1 && ($filters_old[0]=="" || $filters_old[0]=="." || $filters_old[0]=="PASS")) $filters_old = array();

	//set filter entries
	$filters = array_unique(array_merge($filters, $filters_old));
	if(empty($filters))	$filters[] = "PASS";
	$in_file->set($i, 6, implode(";", $filters));
}

//add filter descriptions to header
$filters = array(
	"lt-3-reads" => "Less than 3 supporting tumor reads (filter_vcf).",
	"depth-tum" => "Sequencing depth tumor is less than 20 (filter_vcf).",
	"depth-nor" => "Sequencing depth normal is less than 20 (filter_vcf).",
	"freq-tum" => "Allele frequencies in tumor < 0.05 (filter_vcf).",
	"freq-nor" => "Allele frequencies in normal > 1/6 of allele frequency in tumor (filter_vcf)."
);
$comments = $in_file->getComments();
foreach($filters as $name => $desc)
{
	$comments[] = "#FILTER=<ID={$name},Description=\"{$desc}\">";
}
$in_file->setComments($comments);

//store output file
$in_file->toTSV($out);