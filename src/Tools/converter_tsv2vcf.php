<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_tsv2vcf", "Generate vcf-files from tsv files.");
$parser->addInfile("in",  "Input file in tsv-format.", false);
$parser->addOutfile("out",  "Output file in vcf-format.", false);
// optional
$parser->addStringArray("info",  "Column names that should be included into the INFO part of the vcf file [space separated].", true);
$parser->addFlag("qci",  "Add allele counts '0,0'.");
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

// header
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
$nl = new Matrix();
$nl->addComment("#fileformat=VCFv4.1");
$nl->addComment("#fileDate=".date("Ymd"));
$nl->addComment("#reference=/mnt/storage2/megSAP/data/genomes/{$build}.fa");
$nl->addComment("#INFO=<ID=variant_id,Number=1,Type=String,Description=\"ID of original variant.\">");
$nl->addComment("#FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
$nl->setHeaders(array("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"));

// variants
$ol = Matrix::fromTSV($in);
$idx_chr = $ol->getColumnIndex("chr");
$idx_start = $ol->getColumnIndex("start");
$idx_end = $ol->getColumnIndex("end");
$idx_ref = $ol->getColumnIndex("ref");
$idx_obs = $ol->getColumnIndex("obs");
$idx_info = array();
if(isset($info))
{
	foreach($info as $i)
	{
		$idx_info[$ol->getColumnIndex($i)] = $i;
	}
}

foreach($idx_info as $idx => $name)
{
	$nl->addComment("#INFO=<ID=$name,Number=1,Type=String,Description=\"UNKNOWN.\">");
}

for($i=0;$i<$ol->rows();++$i)
{	
	$row = $ol->getRow($i);
	
	$chr = $row[$idx_chr];
	$start = $row[$idx_start];
	$end = $row[$idx_end];
	$ref = trim($row[$idx_ref],"-");
	$obs = trim($row[$idx_obs],"-");

	if(strlen($ref)!=1 || strlen($obs)!=1)	// is an indel
	{
		list($chr,$start,$ref,$obs) = indel_for_vcf($build, $chr, $start, $ref, $obs);
	}
	
	// set info field
	$info = "variant_id=".$row[$idx_chr]."_".$row[$idx_start]."_".$row[$idx_end]."_".$row[$idx_ref]."_".$row[$idx_obs];
	foreach($idx_info as $idx => $name)
	{
		$info .= ";".$name."=".$row[$idx];
	}
	
	// add row
	$format = "GT";
	$sample = "1/0";
	if($qci)
	{
		$format .= ":AD";
		$sample .= ":0,0";
	}
	$nl->addRow(array($chr,$start,".",$ref,$obs,".",".",$info,$format,$sample));
}
$nl->toTSV($out);