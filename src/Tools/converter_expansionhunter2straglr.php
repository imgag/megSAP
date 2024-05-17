<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_expansionhunter2straglr", "Converts an ExpansionHunter repeat catalog into the straglr input BED format.");
$parser->addInfile("in",  "Input ExpasionHunter catalog file in JSON format.", false);
$parser->addOutfile("out",  "Output straglr BED file.", false);
extract($parser->parse($argv));

function extract_motive($string)
{
	$motive = explode(")", explode("(", $string)[1])[0];
	//replace ambiguous bases with '*'
	$motive = strtr($motive, array("N"=>"*",  "R"=>"*", "Y"=>"*"));
	return $motive;
}

function get_bed_coordinates($string)
{
	list($chr, $start, $end) = explode("\t", strtr($string, array(":"=>"\t", "-"=>"\t")));
	if (!starts_with($chr, "chr")) $chr = "chr".$chr;
	$start = ((int) $start);
	$end = (int) $end;
	return array($chr, $start, $end);
}

$json_file_content = json_decode(file_get_contents($in),true);

$output = array();
$output[] = "#chr\tstart\tend\trepeat_motive\trepeat_id\trepeat_type\tref_size";

foreach ($json_file_content as $repeat)
{
	$ref_region = $repeat["ReferenceRegion"];
	if (is_array($ref_region))
	{
		//multiple repeats per entry
		//split repeat 
		$repeat_motives = explode("*", strtr($repeat["LocusStructure"], array("+"=>"*")));
		$i = 0;
		foreach ($ref_region as $region)
		{
			$line = get_bed_coordinates($region);
			$line[] = extract_motive($repeat_motives[$i]);
			$line[] = $repeat["VariantId"][$i];
			$line[] = $repeat["VariantType"][$i];
			$line[] = ((int) $line[2] - (int) $line[1]) / strlen(trim($line[3])); //add ref size
			$i++;
			$output[] = implode("\t", $line);
		}
	}
	else
	{
		//single repeat per entry
		$line = get_bed_coordinates($ref_region);
		$line[] = extract_motive($repeat["LocusStructure"]);
		$line[] = $repeat["LocusId"];
		$line[] = $repeat["VariantType"];
		$line[] = ((int) $line[2] - (int) $line[1]) / strlen(trim($line[3])); //add ref size
		$output[] = implode("\t", $line);	
	}
	
}

file_put_contents($out, implode("\n", $output));

//sort output file
$parser->exec(get_path("ngs-bits")."BedSort", "-in {$out} -out {$out}");

?>