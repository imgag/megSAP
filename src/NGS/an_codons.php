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
$parser->addFlag("mark_mutation", "Mark the amino acid exchange.");
$parser->addString("field", "VCF field name.", true, "PEPTIDE");

// extract arguments
extract($parser->parse($argv));

// reference genome
$genome = "/mnt/share/data/genomes/{$build}.fa";

// ensGene table
$fp_ensgene = "/mnt/share/data/dbs/Ensembl/ensGene_GRCh37.tsv.gz";
$ensgene_tsv = Matrix::fromTSV($fp_ensgene);
$ensgene = array_column($ensgene_tsv->getData(), NULL, 1);
// ensGene columns
//[0]: bin
//[1]: name
//[2]: chrom
//[3]: strand
//[4]: txStart
//[5]: txEnd
//[6]: cdsStart
//[7]: cdsEnd
//[8]: exonCount
//[9]: exonStarts
//[10]: exonEnds
//[11]: score
//[12]: name2
//[13]: cdsStartStat
//[14]: cdsEndStat
//[15]: exonFrames

// translation table in plain text
$translation = <<<'TABLE'
FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
---M------**--*----M---------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
TABLE;

// construct mapping triplett -> amino acid (one-letter code)
$translation_map = [];
$lines = explode("\n", $translation);
$width = strlen($lines[0]);
for ($c = 0; $c < $width; ++ $c)
{
	$codon = $lines[2][$c] . $lines[3][$c] . $lines[4][$c];
	$aa = $lines[0][$c];
	$translation_map[$codon] = $aa;
}

// translate nucleotide sequence using translation map
function translate($seq)
{
	global $translation_map;
	
	// split into 3-character parts
	$codons = str_split($seq, 3);
	
	// use map to translate
	$aa = [];
	foreach ($codons as $codon) {
		if (in_array($codon, array_keys($translation_map)))
		{
			$aa[] = $translation_map[$codon];
		}
	}

	return implode("", $aa);
}

// read VCF file
$vcf = Matrix::fromTSV($in);

// iterate over VCF records
for ($r = 0; $r < $vcf->rows(); ++ $r)
{
	$fields = $vcf->getRow($r);
		
	$vcf_chr = $fields[0];
	$vcf_pos = $fields[1];
	$vcf_ref = $fields[3];
	$vcf_alt = $fields[4];
	
	// transcript IDs will be extracted from info field/ANN entry
	$transcripts = [];

	// info field
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
		// split multiple entries in ANN
		$ann = explode(",", $info_map["ANN"]);
		foreach ($ann as $a)
		{
			$parts = explode("|", $a);
			$transcripts[] = $parts[6];
		}
	}

	$add_to_info = [];


	foreach ($transcripts as $transcript) {
		if (! isset($ensgene[$transcript]))
		{
			trigger_error("Transcript {$transcript} not found in ensGene table.", E_USER_NOTICE);
			continue;
		}
		$record = $ensgene[$transcript];

		$chr = $record[2];
		$strand = $record[3];

		$txStart = $record[4] + 1;
		$txEnd = $record[5];

		$cdsStart = $record[6] + 1;
		$cdsEnd = $record[7];
		
		$n_exons = $record[8];

		// comma separated, zero based
		$exonStarts = explode(",", $record[9]);
		$exonEnds = explode(",", $record[10]);
	
		// get unspliced 'transcript sequence', define offset to work on it with
		// original coordinates
		$offset = $txStart;
		list($stdout, $stderr) = $parser->exec(get_path("samtools"),
			"faidx {$genome} {$chr}:{$txStart}-{$txEnd} | tail -n+2 | tr -d '\n'",
				false);

		$seq = $stdout[0];

		if (strlen($seq) == 0)
		{
			var_dump($record);
			trigger_error("Empty sequence returned: {$chr}:{$txStart}-{$txEnd}", E_USER_ERROR);
		}
		
		$seq_mut = substr_replace($seq, $vcf_alt, $vcf_pos - $offset, strlen($vcf_ref));
		
		// extract coding sequence with CDS and exon information
		$cds_arr = [];
		$cdsmut_arr = [];
		
		// position of variant in CDS
		$cds_var_pos = -1;
		
		// offsets for insertions/deletions
		$indel_offset_end = 0;
		$indel_offset_start = 0;
		
		// variant in CDS discovered
		$variant_in_cds = false;
		
		// iterate over exons
		for ($e = 0; $e < $n_exons; ++ $e)
		{
			$e_start = $exonStarts[$e] + 1;
			$e_end = $exonEnds[$e];
			
			// skip if exon is upstream/downstream of CDS
			if ($e_end < $cdsStart || $e_start > $cdsEnd)
			{
				continue;
			}
			
			// splice transcript
			// right-most position of exon-start and cds-start
			// if exon overlaps with cds start, results in cds start, else in exon start
			$start = max($e_start, $cdsStart);
			
			// left-most position of exon-end and cds-end
			$end = min($e_end, $cdsEnd);
			
			// check if variant is in this exon
			// if yes, and variant is a deletion or insertion, add/subtract offset
			$variant_in_exon = $vcf_pos >= $start && $vcf_pos <= $end;
			if ($variant_in_exon)
			{
				// flag that variant is in CDS
				$variant_in_cds = true;
				
				// number of inserted/removed bases
				$indel_offset_end = strlen($vcf_alt) - strlen($vcf_ref);
				
				// if variant in this exon, add the number of bases from exon
				// start to the variant to the CDS-position
				$cds_var_pos += $vcf_pos - $start + 1;
			}
			else if ($end < $vcf_pos)
			{
				// if variant not in this exon and variant is downstream, add
				// this exons length to the CDS-position
				$cds_var_pos += $end - $start + 1;
			}
			
			if ($end - $start >= 0)
			{
				// extract exon sequence from transcript sequence:
				// beginning at exon start, and with exon length
				$exon_seq = substr($seq, $start - $offset, $end - $start + 1);
				$cds_arr[] = $exon_seq;
				
				// take substring from mutated transcript sequence, but
				// honor the offsets
				$cdsmut_arr[] = substr($seq_mut,
					$start - $offset + $indel_offset_start,
					$end - $start + 1  + $indel_offset_end);
				
				// reset end-offset and set start-offset
				if ($indel_offset_end != 0)
				{
					$indel_offset_start = $indel_offset_end;
					$indel_offset_end = 0;
				}
			}
		}
		
		$cds = implode("", $cds_arr);
		$cds_mut = implode("", $cdsmut_arr);
//		$cds_mut_2 = substr_replace($cds, $vcf_alt, $cds_var_pos, strlen($vcf_ref));
		
		// reverse complement for transcripts on minus strand
		if ($strand === "-")
		{
			$cds = rev_comp($cds);
			$cds_mut = rev_comp($cds_mut);
//			$cds_mut_2 = rev_comp($cds_mut_2);
			
			// position in reverse complement
			$cds_var_pos = strlen($cds) - $cds_var_pos + 1;
		}
		
		// amino acid position where variant begins
		// 1-based
		$aa_pos = ceil($cds_var_pos / 3);
		
		// translate CDS
		$pep = translate($cds);
		$pep_mut = translate($cds_mut);
//		$pep_mut_2 = translate($cds_mut_2);
		
//		var_dump($fields);
//		var_dump($record);
//		var_dump($cds);
//		var_dump($cds_mut);
//		var_dump($cds_mut_2);
//		var_dump($pep);
//		var_dump($pep_mut);
//		var_dump($pep_mut_2);
//		var_dump($cds_var_pos);
//		var_dump($aa_pos);
		
		$aa_pos0 = $aa_pos - 1;
		
		// codons up to the specified numbers left of the mutation
		$pre_pos = max(0, $aa_pos0 - $flanking_codons);
		$pre_len = $aa_pos0 - $pre_pos;
		$pre = substr($pep, $pre_pos, $pre_len);
		
		// TODO more than one codon?
		$orig =  substr($pep, $aa_pos0, 1);
		$change =  substr($pep_mut, $aa_pos0, 1);
		
		
		$post_pos = min(strlen($pep) - 1, $aa_pos0 + $flanking_codons);
		$post_len = $post_pos - $aa_pos0;
		$post = substr($pep, $aa_pos0 + 1, $post_len);
		
		if ($variant_in_exon)
		{
			if ($mark_mutation)
			{
				$add_to_info[] = sprintf("%s|%s[%s/%s]%s", $transcript, $pre, $orig, $change, $post);
//				$add_to_info[] = sprintf("%s|%s[%s/%s]%s|{$cds}|{$cds_mut}|{$pep}|{$pep_mut}", $transcript, $pre, $orig, $change, $post);
			}
			else
			{
				$add_to_info[] = sprintf("%s|%s%s%s", $transcript, $pre, $change, $post);
			}
		}
		
		
	}

	if (count($add_to_info) > 0)
	{
		$vcf->set($r, 7, $info . ";{$field}=" . implode(",", $add_to_info));
	}
}

$vcf->addComment("#INFO=<ID={$field},Number=.,Type=String,Description=\"Mutated and surrounding peptides from Ensembl peptide reference\">");
$vcf->toTSV($out);