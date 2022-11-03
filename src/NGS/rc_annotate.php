<?php
/**
 * @page rc_annotate
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("rc_annotate", "Annotates a read count table with attributes from GTF file.");
$parser->addInfile("in", "Input read count file in TSV format", false, true);
$parser->addOutfile("out", "Output TSV file for annotated read counts. Can also be the same file as the input file.", false);
$parser->addString("gtfFile", "GTF file containing feature annotations (for mapping identifiers).", false);
$parser->addString("keyId", "Name of key identifier.", true, "gene_id");
$parser->addString("annotationIds", "Annotation identifiers (comma separated) to use as annotation.", true, "gene_name,gene_biotype");
$parser->addString("column_name", "Column name with values to join on.", true, "gene_id");
$parser->addFlag("ignore_version_suffix", "Ignore the Ensembl version suffix when matching key values.");

extract($parser->parse($argv));

// create mapping from GTF file
$mapping = array();

$handle_gtf = fopen2($gtfFile, "r");

while ($line = fgets($handle_gtf))
{
	// skip empty lines and comments
	if ($line === "" || starts_with($line, "#"))
	{
		continue;
	}

	list($chr,
		$source,
		$feature,
		$start,
		$end,
		$score,
		$strand,
		$frame,
		$attributes) = explode("\t", trim($line));

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

    if (array_key_exists($keyId, $attrs))
    {
        $mapping[$attrs[$keyId]] = $attrs;
    }
}
fclose($handle_gtf);

// map function to map key to specified annotation value
function map_id(&$item, $key, $ann_id)
{
	global $mapping;
    global $ignore_version_suffix;

    if ($ignore_version_suffix)
    {
        $item = preg_replace("/\.[0-9]+$/", "", $item);
    }

	if (isset($mapping[$item]) && isset($mapping[$item][$ann_id]))
	{
		$item = $mapping[$item][$ann_id];
	}
	else
	{
		$item = "";
	}
}

// add column with annotation data
$counts = Matrix::fromTSV($in);
$ids = $counts->getCol($counts->getColumnIndex($column_name));

foreach (explode(",", $annotationIds) as $ann_id)
{
	$annotation_column = $ids;
	array_walk($annotation_column, "map_id", $ann_id);
	$counts->removeColByName($ann_id);
	$counts->addCol($annotation_column, $ann_id);
}

$counts->toTSV($out);

?>