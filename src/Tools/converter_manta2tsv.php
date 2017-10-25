<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_manta2tsv", "Generate table from Manta structural variant vcf file.");
$parser->addInfile("in",  "Input file in Manta variant call format.", false);
$parser->addOutfile("out",  "Output file in tabular format.", true);
$parser->addString("tumor_id", "Tumor processed sample identifier.", false);
extract($parser->parse($argv));

$in_tmp = $parser->tempFile();
$parser->exec("gzip -cd", "{$in} > {$in_tmp}", true);

$vcf = Matrix::fromTSV($in_tmp);
//$vcf->printInfo();

if ($vcf->getHeader(9) == $tumor_id)
{
	$tumor_col = 9;
	$normal_col = 10;
}
else if ($vcf->getHeader(10) == $tumor_id)
{
	$tumor_col = 10;
	$normal_col = 9;
}
else
{
	trigger_error("Tumor processed sample identifier is not contained in input file!", E_USER_ERROR);
}


$variants = array();

//print_r($vcf->getComments());

// parse INFO field, i.e. key=value;key=value;...
function parse_info_field($info)
{
	preg_match_all("/([^;= ]+)=([^;= ]+)/", $info, $r); 
	$result = array_combine($r[1], $r[2]);
	return $result;
}

// parse FORMAT and sample field, i.e. PR:SR
function parse_format_field($format, $values)
{
	$keys = explode(":", $format);
	$values_arr = explode(":", $values);
	$result = array_combine($keys, $values_arr);
	return $result;
}

// parse observation "ref,alt", i.e. 56,2 
function parse_observation($str)
{
	list($ref, $alt) = explode(",", $str);
	return array("ref" => $ref, "alt" => $alt);
}

for ($rowidx = 0; $rowidx < $vcf->rows(); $rowidx++)
{
	$row = $vcf->getRow($rowidx);
	
	$chr = $row[0];
	$start = $row[1];
	$id = $row[2];
	$ref = $row[3];
	$alt = $row[4];
	$qual = $row[4];
	$filter = $row[6];
	$info = parse_info_field($row[7]);
	
	$tumor_values = parse_format_field($row[8], $row[$tumor_col]);
	$normal_values = parse_format_field($row[8], $row[$normal_col]);
	
	$variants[$id] = array(
		"type" => $info["SVTYPE"],
		
		"chr" => $chr,
		"start" => $start,
		"end" => array_key_exists("END", $info) ? $info["END"] : ".",
		"filter" => $filter,
		"length" => array_key_exists("SVLEN", $info) ? $info["SVLEN"] : ".",
		
		"score" => $info["SOMATICSCORE"],
		
		"mate_id" => array_key_exists("MATEID", $info) ? $info["MATEID"] : ".",
		"mate_chr" => ".",
		"mate_pos" => ".",
		"mate_filter" => ".",
		
		"tumor_PR_ref" => "n/a",
		"tumor_PR_alt" => "n/a",
		"tumor_SR_ref" => "n/a",
		"tumor_SR_alt" => "n/a",
		
		"normal_PR_ref" => "n/a",
		"normal_PR_alt" => "n/a",
		"normal_SR_ref" => "n/a",
		"normal_SR_alt" => "n/a",
		
		"tumor_PR_freq" => "n/a",
		"tumor_PR_depth" => "n/a",
		"tumor_SR_freq" => "n/a",
		"tumor_SR_depth" => "n/a",
		
		"normal_PR_freq" => "n/a",
		"normal_PR_depth" => "n/a",
		"normal_SR_freq" => "n/a",
		"normal_SR_depth" => "n/a",
		);
	
	if (array_key_exists("PR", $tumor_values))
	{
		$n = parse_observation($tumor_values["PR"]);
		$variants[$id]["tumor_PR_ref"] = $n["ref"];
		$variants[$id]["tumor_PR_alt"] = $n["alt"];
		
		$variants[$id]["tumor_PR_depth"] = $n["ref"] + $n["alt"];
		$variants[$id]["tumor_PR_freq"] = $n["alt"] / ($n["ref"] + $n["alt"]);
	}
	
	if (array_key_exists("SR", $tumor_values))
	{
		$n = parse_observation($tumor_values["SR"]);
		$variants[$id]["tumor_SR_ref"] = $n["ref"];
		$variants[$id]["tumor_SR_alt"] = $n["alt"];
		
		$variants[$id]["tumor_SR_depth"] = $n["ref"] + $n["alt"];
		$variants[$id]["tumor_SR_freq"] = $n["alt"] / ($n["ref"] + $n["alt"]);
	}
	
	if (array_key_exists("PR", $normal_values))
	{
		$n = parse_observation($normal_values["PR"]);
		$variants[$id]["normal_PR_ref"] = $n["ref"];
		$variants[$id]["normal_PR_alt"] = $n["alt"];
		
		$variants[$id]["normal_PR_depth"] = $n["ref"] + $n["alt"];
		$variants[$id]["normal_PR_freq"] = $variants[$id]["normal_PR_depth"] != 0 ? $n["alt"] / ($n["ref"] + $n["alt"]) : "n/a";
	}
	
	if (array_key_exists("SR", $normal_values))
	{
		$n = parse_observation($normal_values["SR"]);
		$variants[$id]["normal_SR_ref"] = $n["ref"];
		$variants[$id]["normal_SR_alt"] = $n["alt"];
		
		$variants[$id]["normal_SR_depth"] = $n["ref"] + $n["alt"];
		$variants[$id]["normal_SR_freq"] = $variants[$id]["normal_SR_depth"] != 0 ? $n["alt"] / ($n["ref"] + $n["alt"]) : "n/a";
	}

}

$table = new Matrix();
$table->setHeaders(array_keys($variants[$vcf->get(0, 2)]));

$resolved_ids = array();

foreach ($variants as $id => $fields)
{
	if (in_array($id, $resolved_ids))
	{
		continue;
	}
	
	if ($fields["mate_id"] !== ".")
	{
		$fields["mate_chr"] = $variants[$fields["mate_id"]]["chr"];
		$fields["mate_pos"] = $variants[$fields["mate_id"]]["start"];
		$fields["mate_filter"] = $variants[$fields["mate_id"]]["filter"];
		$resolved_ids[] = $fields["mate_id"];
	}
	
	$resolved_ids[] = $id;
	$table->addRow($fields);
}

$table->removeCol($table->getColumnIndex("mate_id"));

$table->removeCol($table->getColumnIndex("tumor_PR_ref"));
$table->removeCol($table->getColumnIndex("tumor_PR_alt"));
$table->removeCol($table->getColumnIndex("tumor_SR_ref"));
$table->removeCol($table->getColumnIndex("tumor_SR_alt"));

$table->removeCol($table->getColumnIndex("normal_PR_ref"));
$table->removeCol($table->getColumnIndex("normal_PR_alt"));
$table->removeCol($table->getColumnIndex("normal_SR_ref"));
$table->removeCol($table->getColumnIndex("normal_SR_alt"));

if (isset($out))
{
	$table->toTSV($out);
}
else
{
	$table->toTSV("php://stdout");
}