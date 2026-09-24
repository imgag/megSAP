<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line_straglr arguments
$parser = new ToolBase("converter_expansionhunter2straglr_trgt", "Converts an ExpansionHunter repeat catalog into the straglr input BED format.");
$parser->addInfile("in",  "Input ExpasionHunter catalog file in JSON format.", false);
$parser->addOutfile("out_straglr",  "Output straglr BED file.", false);
$parser->addOutfile("out_trgt",  "Output trgt BED file.", false);
extract($parser->parse($argv));

function extract_motif($string, $replace_ambigious_bases=true)
{
	$motive = explode(")", explode("(", $string)[1])[0];
	//replace ambiguous bases with '*'
	if ($replace_ambigious_bases) $motive = strtr($motive, array("N"=>"*",  "R"=>"*", "Y"=>"*"));
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

function expandDNA(string $sequence)
{
    $map = [
        'A' => ['A'],
        'C' => ['C'],
        'G' => ['G'],
        'T' => ['T'],
        'R' => ['A', 'G'],
        'Y' => ['C', 'T'],
        'S' => ['G', 'C'],
        'W' => ['A', 'T'],
        'K' => ['G', 'T'],
        'M' => ['A', 'C'],
        'B' => ['C', 'G', 'T'],      // not A
        'D' => ['A', 'G', 'T'],      // not C
        'H' => ['A', 'C', 'T'],      // not G
        'V' => ['A', 'C', 'G'],      // not T
        'N' => ['A', 'C', 'G', 'T'], // any
    ];

    $sequence = strtoupper($sequence);
    $results = [''];

    foreach (str_split($sequence) as $char) 
	{
        if (!isset($map[$char])) trigger_error("Unknown IUPAC code: $char", E_USER_ERROR);

        $next = [];
        foreach ($results as $prefix) 
		{
            foreach ($map[$char] as $base) 
			{
                $next[] = $prefix.$base;
            }
        }
        $results = $next;
    }

    return $results;
}

$json_file_content = json_decode(file_get_contents($in),true);

$output_straglr = [];
$output_straglr[] = "#chr\tstart\tend\trepeat_motif\trepeat_id\trepeat_type\tref_size\tref_motif";
$output_trgt = [];

foreach ($json_file_content as $repeat)
{
	$ref_region = $repeat["ReferenceRegion"];
	if (is_array($ref_region))
	{
		//multiple repeats per entry
		//split repeat 
		$repeat_motifs = explode("*", strtr($repeat["LocusStructure"], array("+"=>"*")));
		$i = 0;
		foreach ($ref_region as $region)
		{
			//straglr
			$line_straglr = get_bed_coordinates($region);
			
			$line_straglr[] = extract_motif($repeat_motifs[$i], true);
			$line_straglr[] = $repeat["VariantId"][$i];
			$line_straglr[] = $repeat["VariantType"][$i];
			$line_straglr[] = ((int) $line_straglr[2] - (int) $line_straglr[1]) / strlen(trim($line_straglr[3])); //add ref size
			$line_straglr[] = extract_motif($repeat_motifs[$i], false); //add unmodified motif (for db lookup)
			$output_straglr[] = implode("\t", $line_straglr);

			//trgt:
			$line_trgt = get_bed_coordinates($region);
			//INFO-like 4th column
			$info_trgt = [];
			$info_trgt[] = "ID=".$repeat["VariantId"][$i];
			$info_trgt[] = "VariantType=".$repeat["VariantType"][$i];
			$info_trgt[] = "RefSize=".((int) $line_straglr[2] - (int) $line_straglr[1]) / strlen(trim($line_straglr[3]));
			$info_trgt[] = "MOTIFS=".implode(",", expandDNA(extract_motif($repeat_motifs[$i], false)));
			$info_trgt[] = "RefMotif=".extract_motif($repeat_motifs[$i], false); //add unmodified motif (for db lookup)
			$info_trgt[] = "STRUC=<TR>";

			$line_trgt[] = implode(";", $info_trgt);
			$output_trgt[] = implode("\t", $line_trgt);

			$i++;
			
		}
	}
	else
	{
		//single repeat per entry
		//straglr
		$line_straglr = get_bed_coordinates($ref_region);
		$line_straglr[] = extract_motif($repeat["LocusStructure"], true);
		$line_straglr[] = $repeat["LocusId"];
		$line_straglr[] = $repeat["VariantType"];
		$line_straglr[] = ((int) $line_straglr[2] - (int) $line_straglr[1]) / strlen(trim($line_straglr[3])); //add ref size
		$line_straglr[] = extract_motif($repeat["LocusStructure"], false); //add unmodified motif (for db lookup)
		$output_straglr[] = implode("\t", $line_straglr);

		//trgt:
		$line_trgt = get_bed_coordinates($ref_region);
		//INFO-like 4th column
		$info_trgt = [];
		$info_trgt[] = "ID=".$repeat["LocusId"];
		$info_trgt[] = "VariantType=".$repeat["VariantType"];
		$info_trgt[] = "RefSize=".((int) $line_straglr[2] - (int) $line_straglr[1]) / strlen(trim($line_straglr[3]));
		$info_trgt[] = "MOTIFS=".implode(",", expandDNA(extract_motif($repeat["LocusStructure"], false)));
		$info_trgt[] = "RefMotif=".extract_motif($repeat["LocusStructure"], false); //add unmodified motif (for db lookup)
		$info_trgt[] = "STRUC=<TR>";

		$line_trgt[] = implode(";", $info_trgt);
		$output_trgt[] = implode("\t", $line_trgt);
	}
	
}

file_put_contents($out_straglr, implode("\n", $output_straglr));
file_put_contents($out_trgt, implode("\n", $output_trgt));

//sort output_straglr file
$out_straglr = realpath2($out_straglr); 
$parser->execApptainer("ngs-bits", "BedSort", "-in {$out_straglr} -out {$out_straglr}", [], [dirname($out_straglr)]);
$out_trgt = realpath2($out_trgt); 
$parser->execApptainer("ngs-bits", "BedSort", "-in {$out_trgt} -out {$out_trgt}", [], [dirname($out_trgt)]);

?>