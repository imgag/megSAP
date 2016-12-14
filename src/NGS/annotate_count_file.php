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
//$parser->addString("mapping_file", "File containing the mapping from transcript identifier to gene name", true, get_path("data_folder")."dbs/UCSC/refGene.txt");
$parser->addString("mapping_file", "File containing the mapping from transcript identifier to gene name", true, get_path("data_folder")."dbs/UCSC/refseq2HGNC_twoColumn.tsv");
$parser->addString("mapping_file_delimiter", "Delimiter used in the mapping file", true, "\t");
$parser->addString("counts_file_delimiter", "Delimiter used in the mapping file", true, "\t");
$parser->addInt("id_col", "Column of the mapping file containing the transcript IDs", true, 1);
$parser->addInt("gene_col", "Column of the mapping file containing the gene IDs", true, 2);

extract($parser->parse($argv));

//extracting sub-directories and generating folder structure
create_path($out);

$parser->log("Starting annotation");

// Associative array containing the mapping
$mapping = array();

// Parse the mapping file and save mapping in the associative array
if (($handle = fopen($mapping_file, "r")) !== FALSE) {
    while (($data = fgetcsv($handle, 0, $mapping_file_delimiter)) !== FALSE) {
        $num = count($data);
		
		// data is 0-indexed, but user input is usually 1-indexed
		$mapping[$data[$id_col-1]] = $data[$gene_col-1];
    }
    fclose($handle);
}                   

// Temporary result file
$tmp_res_file = $parser->tempFile();
touch($tmp_res_file);
$tmp_res_file_handle = fopen($tmp_res_file, "w");
$not_mappable = array();

// Read the counts file line by line, annotate and write to the temporary result file
if (($handle = fopen($in, "r")) !== FALSE) {
    while (($data = fgetcsv($handle, 0, $counts_file_delimiter)) !== FALSE) {
		// Get the transcript ID
		$transcript_id = $data[0];
		
		// Try to map the transcript ID
		if(isset($mapping[$transcript_id])) {
			$gene_id = $mapping[$transcript_id];
			$data[] = $gene_id;
		}
		else {
			$gene_id = "";
			$not_mappable[] = $transcript_id;
		}
	
		// Write the result to the temporary result file
		$str = implode($counts_file_delimiter, $data)."\n";
		fwrite($tmp_res_file_handle, $str);
    }
	
    fclose($handle);
}
fclose($tmp_res_file_handle);

// Move temporary file to correct position
copy2($tmp_res_file, $out);
?>