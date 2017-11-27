<?php

/**
	@page an_codons
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_codons", "Add surrounding codon annotation.");

// mandatory arguments
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);

// optional arguments
$parser->addString("build", "Genome build to use.", true, "GRCh37");
$parser->addInt("flanking_codons", "Number of flanking codons.", true, 5);
$parser->addString("field", "VCF field name.", true, "PEPTIDE");

// extract arguments
extract($parser->parse($argv));

// TODO honor build
// peptide FASTA reference
$ref_pep = "/mnt/share/data/dbs/peptide/Homo_sapiens.GRCh37.75.pep.all.fa";
// ENSP--ENST mapping (from FASTA records)
$pep_to_transcript = "/mnt/share/data/dbs/peptide/Homo_sapiens.GRCh37.75.pep.all.tsv";

// read ENSP to ENST mapping
$map_pep_trans_table = Matrix::fromTSV($pep_to_transcript);
$map_pep_trans = array_column($map_pep_trans_table->getData(), 0, 1);

// create amino acid three letter -- one letter code mapping
$aa_code = [
	"Gly"=>"G",
	"Ala"=>"A",
	"Leu"=>"L",
	"Met"=>"M",
	"Phe"=>"F",
	"Trp"=>"W",
	"Lys"=>"K",
	"Gln"=>"Q",
	"Glu"=>"E",
	"Ser"=>"S",
	"Pro"=>"P",
	"Val"=>"V",
	"Ile"=>"I",
	"Cys"=>"C",
	"Tyr"=>"Y",
	"His"=>"H",
	"Arg"=>"R",
	"Asn"=>"N",
	"Asp"=>"D",
	"Thr"=>"T"
];

// read VCF file
$vcf = Matrix::fromTSV($in);

// iterate over VCF records
for ($r = 0; $r < $vcf->rows(); ++ $r)
{
	$transcripts = [];

	$info = $vcf->get($r, 7);
	$info_splt = explode(";", $info);

	// create map for info records
	$info_map = [];
	foreach ($info_splt as $entry)
	{
		if (strpos($entry, "="))
		{
			list($key, $value) = explode("=", $entry, 2);
		}
		else
		{
			$key = $entry;
			$value = NULL;
		}
		$info_map[$key] = $value;
	}

	// extract nucleotide/amino acid change frmo ANN field
	if (isset($info_map["ANN"]))
	{
		$ann = explode(",", $info_map["ANN"]);
		foreach ($ann as $a)
		{
			$parts = explode("|", $a);
			$transcripts[$parts[6]] = ["HGVS.c" => $parts[9], "HGVS.p" => $parts[10]];

			preg_match('/c\.([0-9]+)(.)>(.)/', $parts[9], $matches);
			$transcripts[$parts[6]]["coding"] = $matches;

			preg_match('/p\.([[:alpha:]]{3})([0-9]+)([[:alpha:]]{3})/', $parts[10], $matches);
			$transcripts[$parts[6]]["protein"] = $matches;
		}
	}

	$add_to_info = [];

	foreach ($transcripts as $transcript => $value)
	{
		if (count($value["protein"]) == 4)
		{

			$pos_aa = $value["protein"][2];
			$start = max(1, $pos_aa - $flanking_codons);
			$end = $pos_aa + $flanking_codons;

			$protein_id = $map_pep_trans[$transcript];

			// get peptide sequence
			$pipeline = [];
			$pipeline[] = [get_path("samtools"), "faidx {$ref_pep} {$protein_id}:{$start}-{$end}"];
			$pipeline[] = ["tail", "-n+2"];
			$pipeline[] = ["tr", "-d '\n'"];
			$st_out = $parser->execPipeline($pipeline, "query-pep-ref");

			// unchanged peptide sequence
			$pep_seq = $st_out[0][0];

			// mutated peptide sequence
			// TODO: currently, 3-letter amino acid is inserted (from ANN)
			$repl_pos = $pos_aa - $start;
			$inserted_aa = $aa_code[$value["protein"][3]];
			$pep_seq_mut = substr_replace($pep_seq, $inserted_aa, $repl_pos, 1);

			// add
			$add_to_info[] = sprintf("%s:%s-%s|%s", $protein_id, $start, $end, $pep_seq_mut);
		}
	}

	if (count($add_to_info) > 0)
	{
		$vcf->set($r, 7, $info . ";{$field}=" . implode(",", $add_to_info));
	}
}

$vcf->addComment("#INFO=<ID={$field},Number=.,Type=String,Description=\"Mutated and surrounding peptides from Ensembl peptide reference\">");
$vcf->toTSV($out);