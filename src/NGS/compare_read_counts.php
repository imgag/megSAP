<?php
/**
 * @page compare_read_counts
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("compare_read_counts", "Compare read counts or fpkm values produced by analyze_rna.");
$parser->addInfile("in1", "Counts file for the tumor sample", false, true);
$parser->addInfile("in2", "Counts file for the reference samples", false, true);
$parser->addOutfile("out", "Output file containing the comparison.", false);

//optional parameters
$parser->addString("method", "The comparison method. Options are: 'fc' or 'log_fc'", true, "fc");
$parser->addString("counts_file_delimiter", "Delimiter used in the mapping file", true, "\t");
$parser->addString("mapping_file", "File containing the mapping from transcript identifier to gene name", true, get_path("data_folder")."dbs/UCSC/refseq2HGNC_twoColumn.tsv");
$parser->addInt("id_col", "Column of the mapping file containing the transcript IDs", true, 1);
$parser->addInt("gene_col", "Column of the mapping file containing the gene IDs", true, 2);

extract($parser->parse($argv));

$parser->log("Starting comparison of counts/fpkm values");

// Read the annotation file
// Parse the mapping file and save mapping in the associative array
if (($handle = fopen($mapping_file, "r")) !== FALSE) {
    while (($data = fgetcsv($handle, 0, "\t")) !== FALSE) {
        $num = count($data);
		
		// data is 0-indexed, but user input is usually 1-indexed
		$mapping[$data[$id_col-1]] = $data[$gene_col-1];
    }
    fclose($handle);
}

// Read the two count/fpkm files
$tumor_data = array();
$normal_data = array();
//ini_set('auto_detect_line_endings', true);
$fh1 = fopen($in1, 'r');
$fh2 = fopen($in2, 'r');
// Throw the header away
fgetcsv($fh1, 0, $counts_file_delimiter);
fgetcsv($fh2, 0, $counts_file_delimiter);
// Read in files
while(($line = fgetcsv($fh1, 0, $counts_file_delimiter)) != false) {
  $tumor_data[$line[0]] = $line[1];
}
while(($line = fgetcsv($fh2, 0, $counts_file_delimiter)) != false) {
  $normal_data[$line[0]] = $line[1];
}
fclose($fh1);
fclose($fh2);

// Generate the results file
$fh = fopen($out, 'w');
// Write header
if($method == "fc") {
	fwrite($fh, 'Geneid\tGeneSymbol\tTumor\tReference\tFold_Change\n');
} else if ($method == "log_fc") {
	fwrite($fh, 'Geneid\tGeneSymbol\tTumor\tReference\tLog2_Fold_Change\n');
} else {
	trigger_error("Unknown comparison method: $method", E_USER_ERROR);
}
// Iterate over the tumor data and look for the corresponding reference data
foreach ($tumor_data as $transcript_id => $tumor_val) {
	$line = array($transcript_id);
	// Try to map the transcript id to a gene symbol
	if(isset($mapping[$transcript_id])) {
		$line[] = $mapping[$transcript_id];
	}
	else {
		$line[] = "";
	}
	// Tumor value
	$line[] = strval($tumor_val); 
	// Normal value
    $normal_val = $normal_data[$transcript_id];
	$line[] = strval($normal_val);
	// Compute the (log) fold change. The fold change is not computed if one of the values is 0.
	if($tumor_val == 0 || $normal_val == 0) {
		$line[] = "NA";
	} else {
		$comparison = $tumor_val / $normal_val;
		if($method == "log_fc") {
			$comparison = log($comparison, 2);
		}
		$line[] = strval($comparison);
	}
	// Write the line to the result file
	$str = implode($counts_file_delimiter, $line)."\n";
	fwrite($fh, $str);
}
fclose($fh);
?>