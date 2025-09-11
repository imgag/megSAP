<?php
/**
 * @page rc_annotate
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("rc_annotate", "Annotates a read count table with attributes from GTF file.");
$parser->addInfile("in", "Input read count file in TSV format", false, true);
$parser->addOutfile("out", "Output TSV file for annotated read counts. Can also be the same file as the input file.", false);
$parser->addString("gtf_file", "GTF file containing feature annotations (for mapping identifiers).", false);
$parser->addString("annotationIds", "Annotation identifiers (comma separated) to use as annotation.", true, "gene_name,gene_biotype"); //TODO remove
$parser->addFlag("ignore_version_suffix", "Ignore the Ensembl version suffix when matching key values.");
extract($parser->parse($argv));

// create mapping from GTF file
$mapping = array();
$handle_gtf = fopen2($gtf_file, "r");
while ($line = fgets($handle_gtf))
{
	// skip empty lines and comments
	if ($line === "" || starts_with($line, "#")) continue;

	list($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attributes) = explode("\t", trim($line));

    $attrs = [];
    // list of key-value pairs (strings)
    foreach (preg_split("/;( )?/", $attributes, NULL, PREG_SPLIT_NO_EMPTY) as $kv)
    {
		$ret = preg_split("/ /", $kv, 2, PREG_SPLIT_NO_EMPTY);
		if ($ret !== FALSE)
		{
			list($k, $v) = $ret;
			$attrs[$k] = preg_replace('/(^"|"$)/', "", $v);
		}
		else
		{
			var_dump($kv);
		}
    }

    if (array_key_exists("gene_id", $attrs))
    {
        $mapping[$attrs["gene_id"]] = $attrs;
    }
}
fclose($handle_gtf);


$annotation_ids = explode(",", $annotationIds);

// add column with annotation data
$counts = Matrix::fromTSV($in);

// if columns already exist remove them and update:
$new_col_values = [];
foreach($annotation_ids as $new_col)
{
	$idx = $counts->getColumnIndex($new_col, false, false);
	if ($idx != 0) $counts->removeCol($idx);
	$new_col_values[$new_col] = [];
}

$gene_id_idx = $counts->getColumnIndex("gene_id");
for($row_idx=0; $row_idx < $counts->rows(); ++$row_idx)
{
	$row = $counts->getRow($row_idx);
	$gene_id = $row[$gene_id_idx];
	if ($ignore_version_suffix) $gene_id = preg_replace("/\.[0-9]+$/", "", $gene_id);
	foreach($annotation_ids as $ann_id)
	{
		$anno = "";
		if (isset($mapping[$gene_id]) && isset($mapping[$gene_id][$ann_id]))
		{
			$anno = $mapping[$gene_id][$ann_id];
		}
		$new_col_values[$ann_id][] = $anno;
	}
}

foreach($annotation_ids as $new_col)
{
	$counts->addCol($new_col_values[$new_col], $new_col);
}
$counts->toTSV($out);

?>