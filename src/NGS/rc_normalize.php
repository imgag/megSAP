<?php
/**
 * @page rc_normalize
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parameters
$parser = new ToolBase("rc_normalize", "Normalizes raw read counts produced by featureCounts.");
$parser->addInfile("in", "Input raw read counts file.", false, true);
$parser->addOutfile("out", "Output file for normalized read counts.", false);

//optional parameters
$parser->addEnum("method", "The normalization method.", true, array("rpkm", "cpm"), "rpkm");

extract($parser->parse($argv));

//counts per million mapped reads
function cpm($gene_count, $total_reads) {
	$scaling_permillion = $total_reads / 1e6;
	$rpm = $gene_count / $scaling_permillion;
	return $rpm;
}

//reads per kilobase per million mapped reads
function rpkm($gene_length, $gene_count, $total_reads) {
	$rpm = cpm($gene_count, $total_reads);
	$rpkm = $rpm / ($gene_length / 1e3);
	return $rpkm;
}

//total number of assigned reads
//sum of column 7 (read counts for sample)
$total_reads = $parser->exec("tail", "-n+3 $in | awk -F '\t' '{ x = x + $7 } END { print x }'", false);
$total_reads = intval($total_reads[0][0]);

$in_handle = fopen($in, "r");
if ($in_handle===FALSE) trigger_error("Could not open file '$in' for reading!", E_USER_ERROR);
$out_handle = fopen($out, "w");
if ($out_handle===FALSE) trigger_error("Could not open file '$out' for writing!", E_USER_ERROR);

//skip first two lines
fgetcsv($in_handle, 0, "\t");
fgetcsv($in_handle, 0, "\t");
//process lines
while (($data = fgetcsv($in_handle, 0, "\t")) !== FALSE) {
	//identifier from column 1
	$gene_id = $data[0];
	//gene length from column 6
	$gene_length = $data[5];
	//read count for gene from column 7
	$gene_readcount = $data[6];
	
	//normalized read count
	if ($method == "rpkm") $normalized = rpkm($gene_length, $gene_readcount, $total_reads);
	if ($method == "cpm") $normalized = cpm($gene_readcount, $total_reads);
	
	fwrite($out_handle, $gene_id."\t".$normalized."\n");
}

fclose($in_handle);
fclose($out_handle);

?>