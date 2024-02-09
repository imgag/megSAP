<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("db_converter_ncg", "Converts NCG cancerdrivers_annotation_supporting_evidence.tsv database file (TSV) two seperated tabix-indexed VCF.GZ.");
// optional
$parser->addInfile("in",  "Input file in TSV-format with COSMIC CMC scores. ('-' for STDIN)", false, true);
$parser->addString("outfolder", "Output for the three output files: tumor suppresor gene list (.TXT), the gene region for the tsg (.BED), the reduced db file (.TSV) for oncogenes and tsgs.", false);
$parser->addString("prefix", "Prefix to allow the files to be versioned.", false);
extract($parser->parse($argv));


//compare onco_lines by gene
function cmp($a, $b)
{
	return strcmp($a[1], $b[1]);
}

$tsg_lines = array();
$onco_lines = array();

$header = "Header not set";

foreach (file($in) as $line)
{
	// list(,$gene,,,,$tsg)
	$parts = explode("\t",$line);
	$entrez = $parts[0];
	$gene = $parts[1];
	$cgc_anno = $parts[9];
	$vogel_anno = $parts[10];
	$NCG_onco = trim($parts[12]);
	$NCG_tsg = trim($parts[13]);
	
	$fields = [$entrez, $gene, $cgc_anno, $vogel_anno, $NCG_onco, $NCG_tsg];
	if ($line[0] == "e") # if downloaded by hand the header line doesn't contain the '#' at the start
	{
		$header = "#".implode("\t", $fields);
		continue;
	}
	if($line[0] == "#")
	{
		$header = implode("\t", $fields);
		continue;
	}
	
	if($NCG_onco == "") continue;
	if (in_array($fields, $onco_lines)) continue;
	$onco_lines[] = $fields;
	
	if($NCG_tsg == 1)
	{
		$tsg_lines[] = $gene;
	}
}
//sort onco lines by gene:
usort($onco_lines, "cmp");
$onco_lines_merged = array();

foreach ($onco_lines as $line_parts)
{
	$onco_lines_merged[] = implode("\t", $line_parts);
}


//sort + uniq for tsg:
sort($tsg_lines);
$tsg_lines = array_unique($tsg_lines);

$tsg_file = $outfolder.$prefix."_tsg.txt";
file_put_contents($tsg_file, implode("\n", $tsg_lines));
file_put_contents($outfolder.$prefix."_oncogene.tsv", $header."\n".implode("\n", $onco_lines_merged));

//generate tsg BED file
$pipeline = array();
$pipeline[] = array(get_path("ngs-bits")."GenesToBed", "-source ensembl -mode exon -in $tsg_file");
$pipeline[] = array(get_path("ngs-bits")."BedSort", "-uniq");
$pipeline[] = array(get_path("ngs-bits")."BedMerge", "-merge_names -out ".$outfolder."somatic_tmb_tsg.bed");
$parser->execPipeline($pipeline, "Generate tsg bed file.");



?>
