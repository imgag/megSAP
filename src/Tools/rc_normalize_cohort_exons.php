<?php
/**
 * @page rc_normalize
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parameters
$parser = new ToolBase("rc_normalize_cohort_exons", "Convert and normalize raw exon read counts produced by featureCounts into TSV format.");
$parser->addString("ps_name", "Processed sample name. Needs to be a header of cohort.", true, true);
$parser->addInfile("in_cohort", "Input batch corrected raw gene read counts file.", false, true);
$parser->addInfile("uncorrected_counts", "Input file for uncorrected exon read counts.", false);
$parser->addInfile("normalized_genes", "Raw gene counts for the sample. Length and raw count are necessary for normalization using SRPB.", false);
$parser->addOutfile("out_cohort", "Output file for the normalized read counts for the whole cohort.", false);
$parser->addOutfile("out_sample", "Output file for the normalized read counts for the sample.", false);

//optional parameters
$parser->addString("method_cohort", "The normalization method(s) of exon-level counts, comma-separated if multiple.", true, "srpb");
$parser->addString("method_sample", "The normalization method(s) of exon-level counts, comma-separated if multiple.", true, "raw,rpb,srpb");

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

function getExonColumns($file)
{
	$raw_exon = Matrix::fromTSV($file);
	$raw_exon->removeRow(0);
	$raw_exon->setHeaders([ "gene_id", "chr", "start", "end", "strand", "length", "count" ]);

	// exon lengths
	$lengths = $raw_exon->getCol($raw_exon->getColumnIndex("length"));

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
	
	return [$ids, $ids_gene, $ids_exon, $lengths];
}

function getSRPBscalingFaktor($normalized_genes)
{
	$nrom_genes = Matrix::fromTSV($normalized_genes);
	
	$gene_counts = $nrom_genes->getCol($nrom_genes->getColumnIndex("raw"));
	$gene_lengths = $nrom_genes->getCol($nrom_genes->getColumnIndex("length"));
	
	$rpk_values = [];
	for($idx=0; $idx<count($gene_counts); ++$idx)
	{
		$rpk_values[] = $gene_counts[$idx] / $gene_lengths[$idx] * 1000;
	}
	
	
	return (array_sum($rpk_values) / 1e5);
	
}

//EXON normalization
function rpb($count, $length)
{
	$rpb = [];
	for ($idx = 0; $idx < count($count); ++ $idx)
	{
		$rpb[] = $count[$idx] / $length[$idx];
	}
	
	return $rpb;
}

function srpb($rpb, $scaling_srpb)
{	
	$srpb = [];
	foreach ($rpb as $value)
	{
		$srpb[] = $value / $scaling_srpb;
	}
	return $srpb;
}

function normalize_column($methods, $out_table, $counts, $lengths, $scaling_srpb, $prefix="")
{	
	if (in_array("raw", $methods))
	{
		$out_table->addCol($counts, $prefix."raw", "Number of reads overlapping exon");
	}
	
	if (in_array("rpb", $methods) || in_array("srpb", $methods))
	{
		$rpb = rpb($counts, $lengths);
		
		if (in_array("rpb", $methods))
		{
			$out_table->addCol($rpb, $prefix."rpb", "Reads overlapping exon per base");
		}
		
		if (in_array("srpb", $methods))
		{
			$out_table->addCol(srpb($rpb, $scaling_srpb), $prefix."srpb");
		}
	}
}


$methods_cohort = parse_methods($method_cohort);
$methods_sample = parse_methods($method_sample);

#scaling faktor for srpb: SUM(GENE_COUNTS) / (SUM(GENE_LENGTHS)/100e6)
$scaling_srpb = 0;
if (in_array("srpb", $methods_cohort) || in_array("srpb", $methods_sample)) $scaling_srpb = getSRPBscalingFaktor($normalized_genes);
$parser->log("SRPB scaling factor: $scaling_srpb\n");
 
//output table
$out_sample_table = new Matrix();
$out_cohort_table = new Matrix();

list($ids, $gene_ids, $exon_ids, $exon_lengths) = getExonColumns($uncorrected_counts);

$out_sample_table->addCol($ids, "gene_exon", "Gene exon identifier");
$out_sample_table->addCol($gene_ids, "gene_id", "Gene identifier");
$out_sample_table->addCol($exon_ids, "exon", "Exon coordinates");
$out_sample_table->addCol($exon_lengths, "length", "Exon length");

$out_cohort_table->addCol($ids, "gene_exon", "Gene exon identifier");
$out_cohort_table->addCol($exon_lengths, "length", "Exon length");

$corrected_cohort_counts = Matrix::fromTSV($in_cohort);

for($sample_col=1; $sample_col < $corrected_cohort_counts->cols(); ++$sample_col)
{
	$sample = $corrected_cohort_counts->getHeader($sample_col);
	$sample_counts = $corrected_cohort_counts->getCol($sample_col);
	
	normalize_column($methods_cohort, $out_cohort_table, $sample_counts, $exon_lengths, $scaling_srpb, $sample."_");
	
	if ($sample == $ps_name) normalize_column($methods_sample, $out_sample_table, $sample_counts, $exon_lengths, $scaling_srpb);
}

$out_sample_table->toTSV($out_sample);

//compress output
$tmp_cohort_file = $parser->tempFile("_rc_normalize_cohort_out.tsv");
$out_cohort_table->toTSV($tmp_cohort_file);
$parser->execApptainer("htslib", "bgzip", "-c $tmp_cohort_file > $out_cohort", [], [dirname($out_cohort)]);

?>