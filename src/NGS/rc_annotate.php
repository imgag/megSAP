<?php
/**
 * @page rc_annotate
 * 
 * TODO:
 * - support multiple annotation values
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("rc_annotate", "Annotates a read count file with gene names.");
$parser->addInfile("in", "Input read count file in TSV format", false, true);
$parser->addOutfile("out", "Output TSV file for annotated read counts. Can also be the same file as the input file.", false);
$parser->addString("gtfFile", "GTF file containing feature annotations (for mapping identifiers).", true, get_path("data_folder")."dbs/gene_annotations/GRCh37.gtf");
$parser->addString("keyId", "Name of key identifier.", true, "gene_id");
$parser->addString("annotationId", "Identifier of annotation to add.", true, "gene_name");

$parser->addFlag("ignore_version_suffix", "Ignore the Ensembl version suffix when matching key values.");
$parser->addString("column_name", "Column name with values to join on.", true, "gene_id");

extract($parser->parse($argv));

//read keyId (gene_id) -> annotationId (gene_name) mapping from GTF file
$mapping = array();
$handle_gtf = fopen($gtfFile, "r");
if ($handle_gtf === FALSE)
{
	trigger_error("Could not open file '$gtfFile' for reading!", E_USER_ERROR);
}
while (!feof($handle_gtf))
{
	$line = trim(fgets($handle_gtf));
	if ($line === "")
	{
		continue;
	}

	$parts = explode("\t", $line);

	//parse last column (annotations)
	$annotations = array();
	if (isset($parts[8]))
	{
		foreach (explode(";", $parts[8]) as $anno)
		{
			$anno = trim($anno);
			if (!isset($anno) || $anno === "")
			{
				continue;
			}

			list($key, $value) = explode(" ", $anno);
			$value = trim(substr($value, 1, -1));
			if (isset($key) && isset($value))
			{
				$annotations[$key] = $value;
			}
		}

		if (array_key_exists($keyId, $annotations) && array_key_exists($annotationId, $annotations))
		{
			$mapping[$annotations[$keyId]] = $annotations[$annotationId];
		}
	}
}
fclose($handle_gtf);

// function to query the mapping
function map_annotate($id) {
	global $mapping;
	global $ignore_version_suffix;
	
	if ($ignore_version_suffix)
	{
		$id_nosuffix = preg_replace("/\.[0-9]+$/", "", $id);
		return(isset($mapping[$id_nosuffix]) ? $mapping[$id_nosuffix] : "");
	}
	return(isset($mapping[$id]) ? $mapping[$id] : "");
}

// add column with annotation data
$counts = Matrix::fromTSV($in);
$gene_ids = $counts->getCol($counts->getColumnIndex($column_name));
$annotation_column = array_map("map_annotate", $gene_ids);
$counts->addCol($annotation_column, $annotationId);
$counts->toTSV($out);
