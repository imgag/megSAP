<?php
/**
 * @page rc_normalize
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parameters
$parser = new ToolBase("rc_normalize", "Convert and normalize raw read counts produced by featureCounts into TSV format.");
$parser->addInfile("in", "Input raw read counts file.", false, true);
$parser->addOutfile("out", "Output file for normalized read counts.", false);

//optional parameters
$parser->addString("method", "The normalization method(s), comma-separated if multiple.", true, "raw,cpm,fpkm,tpm");

extract($parser->parse($argv));

$methods = explode(",", $method);

// read raw input into matrix
$raw = Matrix::fromTSV($in);
// remove columns 2,3,4,5
$raw->removeCol(1);
$raw->removeCol(1);
$raw->removeCol(1);
$raw->removeCol(1);
// remove featureCounts header
$raw->removeRow(0);
$raw->setHeaders([ "gene_id", "gene_length", "raw" ]);

// raw counts
$col_raw = $raw->getCol($raw->getColumnIndex("raw"));

// gene length
$col_gene_length = $raw->getCol($raw->getColumnIndex("gene_length"));

// cpm calculation (needed for cpm and fpkm)
if (in_array("cpm", $methods) || in_array("fpkm", $methods))
{
	// scale raw counts with total raw counts in million
	$cpm = [];
	$scaling = array_sum($col_raw) / 1e6;
	foreach ($col_raw as $value)
	{
		$cpm[] = $value / $scaling;
	}
}

// cpm
if (in_array("cpm", $methods))
{
	$raw->addCol($cpm, "cpm");
}

// fpkm
if (in_array("fpkm", $methods))
{
	// scale cpm with gene length in kb
	$fpkm = [];
	for ($idx = 0; $idx < $raw->rows(); ++ $idx)
	{
		$fpkm[] = $cpm[$idx] / ($col_gene_length[$idx] / 1e3);
	}
		
	$raw->addCol($fpkm, "fpkm");
}

// tpm
if (in_array("tpm", $methods))
{
	// scale raw counts with gene length in kb
	$rpk = [];
	for ($idx = 0; $idx < $raw->rows(); ++ $idx)
	{
		$rpk[] = $col_raw[$idx] / ($col_gene_length[$idx] / 1e3);
	}
		
	$tpm = [];
	$scaling = array_sum($rpk) / 1e6;
	foreach ($rpk as $value)
	{
		$tpm[] = $value / $scaling;
	}
	$raw->addCol($tpm, "tpm");
}

// raw
if (! in_array("raw", $methods))
{
	$raw->removeCol($raw->getColumnIndex("raw"));
}

$raw->removeCol($raw->getColumnIndex("gene_length"));
$raw->setComments([]);
$raw->toTSV($out);