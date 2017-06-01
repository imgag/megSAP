<?php
/**
 * @page annotate_count_file
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("annotate_count_file", "Annotates a read count file with gene names.");
$parser->addInfile("in", "Input read count file in TSV format", false, true);
$parser->addOutfile("out", "Output TSV file for annotated read counts. Can also be the same file as the input file.", false);
$parser->addString("gtfFile", "GTF file containing feature annotations (for mapping identifiers).", true, get_path("data_folder")."dbs/gene_annotations/GRCh37.gtf");
$parser->addString("keyId", "Name of key identifier.", true, "gene_id");
$parser->addString("annotationId", "Identifier of annotation to add.", true, "gene_name");
extract($parser->parse($argv));

//read gene_id>gene_name mapping from GTF file
$mapping = array();
$handle_gtf = fopen($gtfFile, "r");
if ($handle_gtf === FALSE) trigger_error("Could not open file '$gtfFile' for reading!", E_USER_ERROR);
while(!feof($handle_gtf))
{
	$line = trim(fgets($handle_gtf));
	if ($line=="") continue;
	
	$parts = explode("\t", $line);
	
	//parse last column (annotations)
	$annotations = array();
	if (isset($parts[8])) {
	foreach(explode(";", $parts[8]) as $anno)
	{
		$anno = trim($anno);
		if (!isset($anno) || $anno=="") continue;
		
		list($key, $value) = explode(" ", $anno);
		$value = trim(substr($value, 1, -1));
		if (isset($key) && isset($value)) {
			$annotations[$key] = $value;
		}
	}
	if (array_key_exists($keyId, $annotations) && array_key_exists($annotationId, $annotations)) {
		$mapping[$annotations[$keyId]] = $annotations[$annotationId];
	}
	}
}
fclose($handle_gtf);

// Read the counts file line by line, annotate and write to the temporary result file
$tmp_res_file = $parser->tempFile();
$tmp_res_file_handle = fopen($tmp_res_file, "w");
$handle = fopen($in, "r");
if ($handle === FALSE) trigger_error("Could not open file '$in' for reading!", E_USER_ERROR);
while ($data = fgetcsv($handle, 0, "\t"))
{
	$transcript_id = $data[0];
	$data[] = isset($mapping[$transcript_id]) ? $mapping[$transcript_id] : "";
	fwrite($tmp_res_file_handle, implode("\t", $data)."\n");
}
fclose($handle);
fclose($tmp_res_file_handle);

// Move temporary file to correct position
copy2($tmp_res_file, $out);
?>