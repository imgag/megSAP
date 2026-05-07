<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_cancerhotspots_check_aminoacid", "");
$parser->addInfile("hotspot_table", "table with genes and residues to check", false);
$parser->addInfile("sequences", "table with genes and full aa sequence", false);
$parser->addString("out", "writes the necessary columns for the NGSD DB if given", true, "");
extract($parser->parse($argv));


$out_lines = [];

$protein_sequences = [];
$geneToTranscript= [];
foreach(file($sequences) as $line)
{
	if (trim($line) == "" || $line[0] == "#") continue;
	
	list($gene, $transcript, $aas) = explode("\t", trim($line));
	
	$protein_sequences[$gene] = $aas;
	$geneToTranscript[$gene] = $transcript;
}

foreach(file($hotspot_table) as $line)
{
	if (trim($line) == "" || $line[0] == "#") continue;
	
	$parts = explode("\t", trim($line));
	
	$gene = $parts[0];
	$residue = $parts[1];
	$variants = $parts[3];
	$samples = $parts[5];
	
	$ref_aa = $residue[0];
	$aa_pos = intval(substr($residue, 1));
	
	if ($protein_sequences[$gene][$aa_pos-1] != $ref_aa)
	{
		trigger_error("Reference amino acid does not match the amino acid in the given sequence for the gene: $gene - $residue", E_USER_ERROR);
	}
		
	if ($out != "")
	{
		foreach(explode("|", $variants) as $var)
		{
			list($alt, $alt_count) = explode(":", $var);
			
			$out_lines[] = implode("\t", [$gene, $geneToTranscript[$gene], $aa_pos, $ref_aa, $alt, $samples, $alt_count]);
		}
		
	}
}

if ($out != "") file_put_contents($out, implode("\n", $out_lines));




?>