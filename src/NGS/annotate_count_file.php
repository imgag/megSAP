<?php
/**
 * @page annotate_count_file
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("annotate_count_file", "Annotates a read count file with gene names.");
$parser->addInfile("in", "Input read count file in tsv format", false, true);
$parser->addOutfile("out", "Output file for annotated read counts. Can also be the same file as the input file.", false);
//optional parameters
$parser->addString("mapping_file", "File containing the mapping from transcript identifier to gene name", true, get_path("data_folder")."dbs/UCSC/refseq2HGNC.tsv");
extract($parser->parse($argv));

//read transcript>gene mapping
$mapping = array();
$handle = fopen($mapping_file, "r");
if ($handle===FALSE) trigger_error("Could not open file '$mapping_file' for reading!", E_USER_ERROR);
while ($data = fgetcsv($handle, 0, "\t"))
{
	$mapping[$data[0]] = $data[1];
}
fclose($handle);                 

// Read the counts file line by line, annotate and write to the temporary result file
$tmp_res_file = $parser->tempFile();
$tmp_res_file_handle = fopen($tmp_res_file, "w");
$handle = fopen($in, "r");
if ($handle===FALSE) trigger_error("Could not open file '$mapping_file' for reading!", E_USER_ERROR);
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