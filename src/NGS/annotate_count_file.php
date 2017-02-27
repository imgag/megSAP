<?php
/**
 * @page annotate_count_file
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("annotate_count_file", "Annotates a read count file with gene names.");
$parser->addInfile("in", "Input read count file in TSV format", false, true);
$parser->addOutfile("out", "Output TSV file for annotated read counts. Can also be the same file as the input file.", false);
extract($parser->parse($argv));

//read gene_id>gene_name mapping from GTF file
$mapping = array();
$mapping_file = get_path("data_folder")."dbs/UCSC/refGene.gtf";
$handle = fopen($mapping_file, "r");
if ($handle===FALSE) trigger_error("Could not open file '$mapping_file' for reading!", E_USER_ERROR);
while(!feof($handle))
{
	$line = trim(fgets($handle));
	if ($line=="") continue;
	
	$parts = explode("\t", $line); 
	
	//parse last column (annotations)
	$annos = array();
	foreach(explode(";", $parts[8]) as $anno)
	{
		$anno = trim($anno);
		if ($anno=="") continue;
		
		list($key, $value) = explode(" ", $anno);
		$annos[$key] = trim(substr($value, 1, -1));
	}
	$mapping[$annos['gene_id']] = $annos['gene_name'];
}
fclose($handle);              

// Read the counts file line by line, annotate and write to the temporary result file
$tmp_res_file = $parser->tempFile();
$tmp_res_file_handle = fopen($tmp_res_file, "w");
$handle = fopen($in, "r");
if ($handle===FALSE) trigger_error("Could not open file '$in' for reading!", E_USER_ERROR);
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