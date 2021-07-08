<?php
/**
 * @page rc_normalize
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parameters
$parser = new ToolBase("rc_normalize", "Convert and normalize raw gene read counts and exon read counts produced by featureCounts into TSV format.");
$parser->addInfile("in", "Input raw gene read counts file.", false, true);
$parser->addInfile("in_exon", "Input raw exon read counts file.", true, true);
$parser->addOutfile("out", "Output file for normalized gene read counts.", false);
$parser->addOutfile("out_exon", "Output file for normalized exon read counts.", true);

//optional parameters
$parser->addString("method", "The normalization method(s), comma-separated if multiple.", true, "raw,cpm,fpkm,tpm");
$parser->addString("method_exon", "The normalization method(s) of exon-level counts, comma-separated if multiple.", true, "raw,rpb,srpb");

extract($parser->parse($argv));

$methods = explode(",", $method);

//output table
$tbl_gene = new Matrix();

// read raw input into matrix
$raw_gene = Matrix::fromTSV($in);
$raw_gene->removeRow(0);
$raw_gene->setHeaders([ "gene_id", "chr", "start", "end", "strand", "length", "count" ]);

// gene IDs
$ids = $raw_gene->getCol($raw_gene->getColumnIndex("gene_id"));
$tbl_gene->addCol($ids, "gene_id");

// gene length
$length = $raw_gene->getCol($raw_gene->getColumnIndex("length"));

// raw counts
$count = $raw_gene->getCol($raw_gene->getColumnIndex("count"));

// per-million scaling
$scaling = array_sum($count) / 1e6;

// cpm calculation (needed for cpm and fpkm)
if (in_array("cpm", $methods) || in_array("fpkm", $methods))
{
	// scale raw counts with total raw counts in million
	$cpm = [];
	foreach ($count as $value)
	{
		$cpm[] = $value / $scaling;
	}
}


// raw
if (in_array("raw", $methods))
{
	$tbl_gene->addCol($count, "raw");
}

// cpm
if (in_array("cpm", $methods))
{
	$tbl_gene->addCol($cpm, "cpm");
}

// fpkm
if (in_array("fpkm", $methods))
{
	// scale cpm with gene length in kb
	$fpkm = [];
	for ($idx = 0; $idx < $raw_gene->rows(); ++ $idx)
	{
		$fpkm[] = $cpm[$idx] / ($length[$idx] / 1e3);
	}
		
	$tbl_gene->addCol($fpkm, "fpkm");
}

// tpm
if (in_array("tpm", $methods))
{
	// scale raw counts with gene length in kb
	$rpk = [];
	for ($idx = 0; $idx < $raw_gene->rows(); ++ $idx)
	{
		$rpk[] = $count[$idx] / ($length[$idx] / 1e3);
	}
		
	$tpm = [];
	$scaling = array_sum($rpk) / 1e6;
	foreach ($rpk as $value)
	{
		$tpm[] = $value / $scaling;
	}
	$tbl_gene->addCol($tpm, "tpm");
}

// write gene-level table
$tbl_gene->toTSV($out);

if (isset($in_exon))
{
	$methods_exon = explode(",", $method_exon);

	$tbl_exon = new Matrix();
	$raw_exon = Matrix::fromTSV($in_exon);
	$raw_exon->removeRow(0);
	$raw_exon->setHeaders([ "gene_id", "chr", "start", "end", "strand", "length", "count" ]);

	// raw read counts per exon
	$count_exon = $raw_exon->getCol($raw_exon->getColumnIndex("count"));
	// exon length
	$length_exon = $raw_exon->getCol($raw_exon->getColumnIndex("length"));

	// use gene_id:chr:exon start - exon end:strand as identifier
	$ids = [];
	$ids_gene = [];
	$ids_exon = [];
	for ($idx = 0; $idx < $raw_exon->rows(); ++ $idx)
	{
		$gene = $raw_exon->get($idx, 0);
		$chr = $raw_exon->get($idx, 1);
		$start = $raw_exon->get($idx, 2);
		$end = $raw_exon->get($idx, 3);
		$strand = $raw_exon->get($idx, 4);
		$ids[] = "{$gene}:{$chr}:{$start}-{$end}:{$strand}";
		$ids_gene[] = $gene;
		$ids_exon[] = "{$chr}:{$start}-{$end}:{$strand}";
	}
	$tbl_exon->addCol($ids, "gene_exon", "Gene exon identifier");
	$tbl_exon->addCol($ids_gene, "gene_id", "Gene identifier");
	$tbl_exon->addCol($ids_exon, "exon", "Exon coordinates");

	if (in_array("raw", $methods_exon))
	{
		$tbl_exon->addCol($count_exon, "raw", "Number of reads overlapping exon");
	}

	if (in_array("rpb", $methods_exon) || in_array("srpb", $methods_exon))
	{
		// scale raw counts with gene length in b
		$rpb = [];
		for ($idx = 0; $idx < $raw_exon->rows(); ++ $idx)
		{
			$rpb[] = $count_exon[$idx] / $length_exon[$idx];
		}
	}

	if (in_array("rpb", $methods_exon))
	{
		$tbl_exon->addCol($rpb, "rpb", "Reads overlapping exon per base");
	}

	if (in_array("srpb", $methods_exon))
	{
		$srpb = [];
		foreach ($rpb as $value)
		{
			$srpb[] = $value / $scaling * 100;
		}

		$tbl_exon->addCol($srpb, "srpb", "Reads overlapping exon per base, scaled by total number of reads in 100 million");
	}

	$tbl_exon->addComment("Total library size (from gene-level counting): " . array_sum($count));
	// write exon-level table
	$tbl_exon->toTSV($out_exon);

}