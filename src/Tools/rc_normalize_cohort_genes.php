<?php
/**
 * @page rc_normalize
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parameters
$parser = new ToolBase("rc_normalize_cohort_genes", "Convert and normalize raw gene read counts produced by featureCounts into TSV format.");
$parser->addString("ps_name", "Processed sample name. Needs to be a header of cohort.", true, true);
$parser->addInfile("in_cohort", "Input raw gene read counts file.", false, true);
$parser->addInfile("uncorrected_counts", "Output file for normalized gene read counts.", false);
$parser->addOutfile("out_cohort", "Output file for the normalized read counts for the whole cohort.", false);
$parser->addOutfile("out_sample", "Output file for the normalized read counts for the sample.", false);

//optional parameters
$parser->addString("method_cohort", "The normalization method(s), comma-separated if multiple.", true, "tpm");
$parser->addString("method_sample", "The normalization method(s), comma-separated if multiple.", true, "raw,cpm,fpkm,tpm");

extract($parser->parse($argv));

function parse_methods($method_arg)
{
	if (contains($method_arg, ","))
	{
		return explode(",", $method_arg);
	} 
	else
	{
		return [$method_arg];
	}
}


function getGeneIDandLength($file)
{
	// read uncorrected input into matrix for gene_id and length
	$raw_gene = Matrix::fromTSV($file);
	$raw_gene->removeRow(0);
	$raw_gene->setHeaders([ "gene_id", "chr", "start", "end", "strand", "length", "count" ]);
	// gene IDs
	$ids = $raw_gene->getCol($raw_gene->getColumnIndex("gene_id"));
	// gene length
	$gene_lengths = $raw_gene->getCol($raw_gene->getColumnIndex("length"));
	
	return [$ids, $gene_lengths];
}


//GENE normalization
function cpm($count)
{
	//scale raw counts with total raw counts in million
	$scaling = array_sum($count) / 1e6;
	$cpm = [];
	foreach ($count as $value)
	{
		$cpm[] = $value / $scaling;
	}
	return $cpm;
}

function fpkm($cpm, $length)
{
	$fpkm = [];
	for ($idx = 0; $idx < count($cpm); ++ $idx)
	{
		$fpkm[] = $cpm[$idx] / ($length[$idx] / 1e3);
	}
		
	return $fpkm;
}

function tpm($count, $length)
{
	$rpk = [];
	for ($idx = 0; $idx < count($count); ++ $idx)
	{
		$rpk[] = $count[$idx] / ($length[$idx] / 1e3);
	}
		
	$tpm = [];
	$rpk_sum = array_sum($rpk);
	foreach ($rpk as $value)
	{
		$tpm[] = $value / $rpk_sum *1e6;
	}
	
	return $tpm;
}

function normalize_column($methods, $out_table, $counts, $lengths, $prefix="")
{
	if (in_array("raw", $methods))
	{
		$out_table->addCol($counts, $prefix."raw");
	}
	
	if (in_array("cpm", $methods) || in_array("fpkm", $methods))
	{
		$cpm = cpm($counts);
		
		if (in_array("cpm", $methods))
		{
			$out_table->addCol($cpm, $prefix."cpm");
		}
		
		if (in_array("fpkm", $methods))
		{
			// scale cpm with gene length in kb
			$out_table->addCol(fpkm($cpm, $lengths), $prefix."fpkm");
		}
	}
	
	if (in_array("tpm", $methods))
	{
		$out_table->addCol(tpm($counts, $lengths), $prefix."tpm");
	}
}


$methods_cohort = parse_methods($method_cohort);
$methods_sample = parse_methods($method_sample);

//output table
$out_sample_table = new Matrix();
$out_cohort_table = new Matrix();

list($ids, $gene_lengths) = getGeneIDandLength($uncorrected_counts);
$out_sample_table->addCol($ids, "gene_id");
$out_sample_table->addCol($gene_lengths, "length");
$out_cohort_table->addCol($ids, "gene_id");
$out_cohort_table->addCol($gene_lengths, "length");

$corrected_cohort_counts = Matrix::fromTSV($in_cohort);

for($sample_col=1; $sample_col < $corrected_cohort_counts->cols(); ++$sample_col)
{
	$sample = $corrected_cohort_counts->getHeader($sample_col);
	$sample_counts = $corrected_cohort_counts->getCol($sample_col);
	
	normalize_column($methods_cohort, $out_cohort_table, $sample_counts, $gene_lengths, $sample."_");
	
	if ($sample == $ps_name) normalize_column($methods_sample, $out_sample_table, $sample_counts, $gene_lengths);
}

$out_sample_table->toTSV($out_sample);
//compress output
$tmp_cohort_file = $parser->tempFile("_rc_normalize_cohort_out.tsv");
$out_cohort_table->toTSV($tmp_cohort_file);
$parser->execApptainer("htslib", "bgzip", "-c $tmp_cohort_file > $out_cohort", [], [dirname($out_cohort)]);

?>