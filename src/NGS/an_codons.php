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

// reference genome
$genome = get_path("local_data")."/{$build}.fa";

// ensGene table
$fp_ensgene = get_path("data_folder")."/dbs/Ensembl/ensGene_GRCh37.tsv.gz";
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
	$fields = $vcf->getRow($r);
		
	$vcf_chr = $fields[0];
	$vcf_pos = $fields[1];
	$vcf_ref = $fields[3];
	$vcf_alt = $fields[4];
	
	// transcript IDs will be extracted from info field/ANN entry
	$transcripts = [];
	$transcripts_snpeff = [];

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
			
			$transcripts_snpeff[$parts[6]] = ["HGVS.c" => $parts[9], "HGVS.p" => $parts[10]];

			preg_match('/c\.([0-9]+)(.)>(.)/', $parts[9], $matches);
			$transcripts_snpeff[$parts[6]]["coding"] = $matches;

			$m = preg_match('/p\.([[:alpha:]]{3})([0-9]+)([[:alpha:]]{3})/', $parts[10], $matches);
			if ($m &&
				in_array($matches[1], array_keys($aa_code)) &&
				in_array($matches[3], array_keys($aa_code)))
			{
				$transcripts_snpeff[$parts[6]]["protein"] = $matches;
			}
			else
			{
				$transcripts_snpeff[$parts[6]]["protein"] = [];
			}
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
		
		$exonFrames = explode(",", preg_replace("/,$/", "", $record[15]));
	
		// skip non-coding transcripts
		if ($cdsStart > $cdsEnd)
		{
			//trigger_error("Skipped transcript {$transcript}, no CDS.", E_USER_NOTICE);
			continue;
		}
		
		// get unspliced 'transcript sequence', define offset to work on it with
		// original coordinates
		$offset = $txStart;
		list($stdout, $stderr) = $parser->exec(get_path("samtools"),
			"faidx {$genome} {$chr}:{$txStart}-{$txEnd} | tail -n+2 | tr -d '\n'",
				false);

		$seq = $stdout[0];

		if (strlen($seq) == 0)
		{
			trigger_error("Empty sequence returned: {$chr}:{$txStart}-{$txEnd}", E_USER_ERROR);
		}
		
		$seq_mut = substr_replace($seq, $vcf_alt, $vcf_pos - $offset, strlen($vcf_ref));
		
		// extract coding sequence with CDS and exon information
		$cds_arr = [];
		$cdsmut_arr = [];
		
		// position of variant in CDS
		$cds_var_pos = 0;
		
		// offsets for insertions/deletions
		$indel_offset_end = 0;
		$indel_offset_start = 0;
		
		// variant in CDS discovered
		$variant_in_cds = false;
		
		// number of exons that are skipped (UTR)
		$skip_exons_5prime = 0;
		$skip_exons_3prime = 0;
		
		// iterate over exons
		for ($e = 0; $e < $n_exons; ++ $e)
		{
			$e_start = $exonStarts[$e] + 1;
			$e_end = $exonEnds[$e];
			
			// skip if exon is upstream/downstream of CDS
			if ($e_end < $cdsStart)
			{
				$skip_exons_5prime += 1;
				continue;
			}
			if ($e_start > $cdsEnd)
			{
				$skip_exons_3prime += 1;
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
		
		if ($strand === "+")
		{
			// first exon frame
			if ($skip_exons_5prime >= count($exonFrames))
			{
				var_dump($transcript);
				var_dump($skip_exons_5prime);
				var_dump($exonFrames);
			}
			$frame = $exonFrames[$skip_exons_5prime];
		}
		
		// reverse complement for transcripts on minus strand
		if ($strand === "-")
		{
			$cds = rev_comp($cds);
			$cds_mut = rev_comp($cds_mut);
			
			// position in reverse complement
			$cds_var_pos = strlen($cds) - $cds_var_pos + 1;
			
			// last exon from table is first exon to be translated 
			$frame = $exonFrames[ $n_exons - 1 - $skip_exons_3prime ];
		}
		
		
		// honor the frame offset
		// frame 2 -> NN<cds>, i.e. remove one nucleotide
		// frame 1 -> N<cds>, i.e. remove two nucleotides
		// frame 0 -> <cds>, as-is
		$frame_offset = (3 - $frame) % 3;
		$cds = substr($cds, $frame_offset);
		$cds_mut = substr($cds_mut, $frame_offset);
		
		// modify cds_var_pos to reflect the same change
		$cds_var_pos -= $frame_offset;
		
		// amino acid position where variant begins
		// 1-based
		$aa_pos = intval(ceil($cds_var_pos / 3));
		
		// translate CDS
		$pep = translate($cds);
		$pep_mut = translate($cds_mut);

		// check nuc pos and aa pos against snpeff
		if (isset($transcripts_snpeff[$transcript]["coding"][1]) &&
			$transcripts_snpeff[$transcript]["coding"][1] != $cds_var_pos)
		{
			trigger_error(sprintf("%s %s %s %s %s variant position in CDS differs: ann=%s / new=%d %s>%s / new+oldpos= %s>%s ; f=%d, f_offset=%d",
				$vcf_chr, $vcf_pos, $vcf_ref, $vcf_alt, $transcript,
				$transcripts_snpeff[$transcript]["coding"][0],
				$cds_var_pos,
				$cds[$cds_var_pos - 1],
				$cds_mut[$cds_var_pos - 1],
				$cds[$transcripts_snpeff[$transcript]["coding"][1] - 1],
				$cds_mut[$transcripts_snpeff[$transcript]["coding"][1] - 1],
				$frame,
				$frame_offset
				),
				E_USER_WARNING);
		}
		if (isset($transcripts_snpeff[$transcript]["protein"][2]) &&
			$transcripts_snpeff[$transcript]["protein"][2] != $aa_pos)
		{
			trigger_error(sprintf("%s %s %s %s %s variant position in peptide differs: ann=%s / new=%d %s > %s / new+oldpos= %s>%s ; f=%d, f_offset=%d",
				$vcf_chr, $vcf_pos, $vcf_ref, $vcf_alt, $transcript,
				$transcripts_snpeff[$transcript]["protein"][0],
				$aa_pos,
				$pep[$aa_pos - 1],
				$pep_mut[$aa_pos - 1],
				$pep[$transcripts_snpeff[$transcript]["protein"][2] - 1],
				$pep_mut[$transcripts_snpeff[$transcript]["protein"][2] - 1],
				$frame,
				$frame_offset
				),
				E_USER_WARNING);
		}
		
		// position of changed amino acid in string (0-based)
		$aa_pos0 = $aa_pos - 1;
		
		// codons up to the specified numbers left of the mutation
		$pre_pos = max(0, $aa_pos0 - $flanking_codons);
		$pre_len = $aa_pos0 - $pre_pos;
		$pre = substr($pep_mut, $pre_pos, $pre_len);
		
		// TODO more than one codon?
		$orig =  substr($pep, $aa_pos0, 1);
		$change =  substr($pep_mut, $aa_pos0, 1);
		
		
		$post_pos = min(strlen($pep_mut) - 1, $aa_pos0 + $flanking_codons);
		$post_len = $post_pos - $aa_pos0;
		$post = substr($pep_mut, $aa_pos0 + 1, $post_len);
		
		if ($variant_in_cds)
		{
			$wt = $pre.$orig.$post;
			$mut = $pre.$change.$post;
			
			if(strpos($mut,"*")!==FALSE)	$mut = substr($mut,0,strpos($mut,"*")+1);
			
			$add_to_info[] = sprintf("%s|%s|%s", $transcript, $wt, $mut);
		}
		else
		{
			if (isset($transcripts_snpeff[$transcript]["protein"][0]))
			{
				var_dump($transcripts_snpeff[$transcript]["protein"]);
				trigger_error(sprintf("Transcript %s not annotated but SnpEff record available!.", $transcript), E_USER_NOTICE);
			}
			
		}
		
		
	}

	if (count($add_to_info) > 0)
	{
		$vcf->set($r, 7, $info . ";{$field}=" . implode(",", $add_to_info));
	}
}

$vcf->addComment("#INFO=<ID={$field},Number=.,Type=String,Description=\"Peptide prediction\">");
$vcf->toTSV($out);