<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_straglr2methylartist", "Converts a straglr repeat expansion catalog into a input file for methylartist");
$parser->addInfile("catalog",  "Straglr catalog in BED format.", false);
$parser->addOutfile("out",  "Output file in TSV-format.", false);
$parser->addInfile("imprinting_tsv", "Additional imprinting TSV file which will be added to the repeat expansions.", true);
$parser->addInfile("promotor_tsv", "Additional promotor TSV file which will be added to the repeat expansions.", true);
extract($parser->parse($argv));


if ($imprinting_tsv != "")
{
	$output = Matrix::fromTSV($imprinting_tsv);

	//basic file check
	if ($output->cols() != 8) trigger_error("Column count of imprinting input file has to be 8 (is: ".$output->cols().")!", E_USER_ERROR);
}
else
{
	$output = new Matrix();
	$output->setHeaders(array("identifier", "title", "gene symbol", "gene chr", "gene start", "gene end", "highlight start", "highlight end"), null);
}

if ($promotor_tsv != "")
{
	$tmp = Matrix::fromTSV($promotor_tsv);

	//basic file check
	if ($tmp->cols() != 8) trigger_error("Column count of promotor input file has to be 8 (is: ".$tmp->cols().")!", E_USER_ERROR);

	//append to list
	for($i = 0; $i < $tmp->rows(); ++$i)
	{
		$output->addRow($tmp->getRow($i));
	}
}



//read repeat expansion catalog
$re_bed = Matrix::fromTSV($catalog);
//basic file check
if ($re_bed->cols() != 8) trigger_error("Column count of straglr catalog file has to be 8 (is: ".$re_bed->cols().")!", E_USER_ERROR);

for($i = 0; $i < $re_bed->rows(); ++$i)
{
	list($re_chr, $re_start, $re_end, $repeat_motif, $repeat_id, $repeat_type, $ref_size, $ref_motif) = $re_bed->getRow($i);
	
	$gene = explode("_", $repeat_id)[0];	

	//debug 
	print $gene."\n";

	//get gene region
	$gene_pipeline = array();
	$gene_pipeline[] = array("echo", $gene);
	$gene_pipeline[] = array("", $parser->execApptainer("ngs-bits", "GenesToBed", "-source ensembl -mode gene", [], [], true));
	$gene_pipeline[] = array("", $parser->execApptainer("ngs-bits", "BedMerge", "", [], [], true));
	list($stdout, $stderr) = $parser->execPipeline($gene_pipeline, "gene region pipeline");
	
	list($gene_chr, $gene_start, $gene_end) = explode("\t", $stdout[0]);

	//sanity check:
	if ($re_chr != $gene_chr) trigger_error("Gene chromosome doesn't match chromosome of repeat ({$gene_chr} vs. {$re_chr})!", E_USER_ERROR);

	$row = array();
	$row[] = "RepeatExpansion_".$repeat_id;
	$row[] = "Repeat Expansion ".$repeat_id;
	$row[] = $gene;
	$row[] = $gene_chr;
	$row[] = $gene_start;
	$row[] = $gene_end;
	$row[] = $re_start;
	$row[] = $re_end;

	$output->addRow($row);
}

$output->toTSV($out);
