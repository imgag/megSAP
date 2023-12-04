<?php

/**
 * @page transcript_variant_primer
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("transcript_variant_primer", "Design primers surrounding a variant in a specific transcript.");
$parser->addString("prefix", "User prefix/name for primer design.", false);
$parser->addString("transcript_id", "Target transcript identifier (Ensembl).", false);
$parser->addInt("variant_position", "Coding (start) position of variant.", false);
$parser->addOutfile("out", "Output table in TSV format.", false);
//optional
$parser->addInt("variant_length", "Variant length.", true, 1);
$parser->addString("build", "Genome build to use.", true, "GRCh38");
$parser->addOutfile("primer3", "Keep primer3 result file.", true);
$parser->addOutfile("bam", "Keep primer alignment BAM file.", true);
$parser->addInt("n", "Number of primer pairs to generate.", true, 3);
$parser->addInt("size_opt", "Optimal primer size.", true, 18);
$parser->addInt("size_min", "Minimum primer size.", true, 15);
$parser->addInt("size_max", "Maximum primer size.", true, 21);
$parser->addInt("product_min", "Minimum product size.", true, 190);
$parser->addInt("product_max", "Maximum product size.", true, 250);
extract($parser->parse($argv));

$reference = "/mnt/storage2/megSAP/data/genomes/${build}_cDNA.fa";
$primer3_config = get_path("primer3")."primer3_config/";

//ensgene table (cDNA, CDS, exons)
$fp_ensgene = get_path("GRCh37_data_folder")."/dbs/Ensembl/ensGene_GRCh37.tsv.gz";
$ensgene_tsv = Matrix::fromTSV($fp_ensgene);
$ensgene = array_column($ensgene_tsv->getData(), NULL, 1);

//offset cDNA/CDS
$plus_strand = $ensgene[$transcript_id][3] === "+";

$exon_starts = explode(",", rtrim($ensgene[$transcript_id][9], ","));
$exon_ends = explode(",", rtrim($ensgene[$transcript_id][10], ","));

$offset = 0;
$n_exons = $ensgene[$transcript_id][8];
$cdsStart = $ensgene[$transcript_id][6] + 1;
$cdsEnd = $ensgene[$transcript_id][7];
// iterate over exons
if ($plus_strand)
{
	for ($e = 0; $e < $n_exons; ++ $e)
	{
		$e_start = $exon_starts[$e] + 1;
		$e_end = $exon_ends[$e];

		// UTR exon - add to offset
		if ($e_start < $cdsStart && $e_end < $cdsStart)
		{
			$offset += $e_end - $e_start + 1;
		}

		// partly UTR
		else if ($e_start < $cdsStart && $e_end >= $cdsStart)
		{
			$offset += $cdsStart - $e_start;
		}
	}
}
else
{
		for ($e = 0; $e < $n_exons; ++ $e)
	{
		$e_start = $exon_starts[$e] + 1;
		$e_end = $exon_ends[$e];

		// UTR exon - add to offset
		if ($e_start > $cdsEnd && $e_end > $cdsEnd)
		{
			$offset += $e_end - $e_start + 1;
		}

		// partly UTR
		else if ($e_end > $cdsEnd && $e_start <= $cdsEnd)
		{
			$offset += $e_end - $cdsEnd;
		}
	}
}

function sub($a, $b)
{
	return $a - $b;
}
if (!$plus_strand)
{
	$exon_ends = array_reverse($exon_ends);
	$exon_starts = array_reverse($exon_starts);
}
$exon_lengths = array_map("sub", $exon_ends, $exon_starts);
$transcript_variant_position = $variant_position + $offset;

//drop last value (end position, i.e. not a junction which can be covered)
$sum = 0;
foreach ($exon_lengths as $l) {
	$sum += $l;
	$junctions_arr[] = $sum;
}
array_pop($junctions_arr);
$junctions = implode(" ", $junctions_arr);

//get transcript (cDNA) sequence
$pipeline = [
	[get_path("samtools"), "faidx $reference $transcript_id"],
	["tail", "-n+2"],
	["tr", "-d '\n'"]
];
list($stdout) = $parser->exec($pipeline, "get transcript cDNA");
$transcript_sequence = $stdout[0];

//primer3
$primer3_input = <<<EOT
SEQUENCE_ID={$transcript_id}
SEQUENCE_TEMPLATE={$transcript_sequence}
SEQUENCE_TARGET={$transcript_variant_position},{$variant_length}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE={$size_opt}
PRIMER_MIN_SIZE={$size_min}
PRIMER_MAX_SIZE={$size_max}
PRIMER_PRODUCT_SIZE_RANGE={$product_min}-{$product_max}
SEQUENCE_OVERLAP_JUNCTION_LIST={$junctions}
PRIMER_NUM_RETURN={$n}
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=0
PRIMER_THERMODYNAMIC_PARAMETERS_PATH={$primer3_config}
=
EOT;
$primer3_in = $parser->tempFile();
file_put_contents($primer3_in, $primer3_input);
$primer3_out = isset($primer3) ? $primer3 : $parser->tempFile();
$parser->exec(get_path("primer3")."primer3_core", "<$primer3_in >$primer3_out", true);

//parse primer3 output into array
$output = [];
$results = [];
foreach (file($primer3_out) as $line)
{
	$line = trim($line);
	if (preg_match("/PRIMER_(PAIR|LEFT|RIGHT)_([0-9]+)_([^=]+)=(.+)/", $line, $matches))
	{
		if (!array_key_exists($matches[2], $results))
		{
			$results[$matches[2]] = [];
		}

		if (!array_key_exists($matches[1], $results[$matches[2]]))
		{
			$results[$matches[2]][$matches[1]] = [];
		}

		$results[$matches[2]][$matches[1]][$matches[3]] = $matches[4];
	}
	else if (preg_match("/PRIMER_(PAIR|LEFT|RIGHT)_([0-9]+)=(.+)/", $line, $matches))
	{
		if (!array_key_exists($matches[2], $results))
		{
			$results[$matches[2]] = [];
		}

		if (!array_key_exists($matches[1], $results[$matches[2]]))
		{
			$results[$matches[2]][$matches[1]] = [];
		}

		$results[$matches[2]][$matches[1]]["position"] = $matches[3];
	}
	else if (preg_match("/([^=]+)=(.+)/", $line, $matches))
	{
		$output[$matches[1]] = $matches[2];
	}
}

//prepare FASTA export of primers
$fasta = [];
foreach ($results as $num => $primer)
{
	$fasta[] = <<<EOT
>primer{$num}_left_{$primer["PAIR"]["PRODUCT_SIZE"]}bp
{$primer["LEFT"]["SEQUENCE"]}
>primer{$num}_right_{$primer["PAIR"]["PRODUCT_SIZE"]}bp
{$primer["RIGHT"]["SEQUENCE"]}
EOT;
}
$primers_fasta = $parser->tempFile(".fa");
file_put_contents($primers_fasta, implode("\n", $fasta));
//$parser->log("Primers FASTA:", file($primers_fasta));

//annotate transcript sequence with primers and variants
$vis = [];
foreach ($results as $num => $primer)
{
	list($pos, $len) = explode(",", $output["SEQUENCE_TARGET"]);
	list($left_pos, $left_len) = explode(",", $primer["LEFT"]["position"]);
	list($right_pos_end, $right_len) = explode(",", $primer["RIGHT"]["position"]);
	// primer3 gives right end of right primer, calculate left position
	$right_pos = $right_pos_end - $right_len + 1;

	$vis[$num] = [];
	// start up to left primer
	$vis[$num][] = substr($output["SEQUENCE_TEMPLATE"], 0, $left_pos);
	// left primer
	$vis[$num][] = substr($output["SEQUENCE_TEMPLATE"], $left_pos, $left_len);
	// after left primer up to variant
	$vis[$num][] = substr($output["SEQUENCE_TEMPLATE"], $left_pos + $left_len, $pos - ($left_pos + $left_len));
	// variant
	$vis[$num][] = substr($output["SEQUENCE_TEMPLATE"], $pos, $len);
	// after variant up to right primer
	$vis[$num][] = substr($output["SEQUENCE_TEMPLATE"], $pos + $len, $right_pos - ($pos + $len));
	// right primer
	$vis[$num][] = substr($output["SEQUENCE_TEMPLATE"], $right_pos, $right_len);
	// after right primer to end
	$vis[$num][] = substr($output["SEQUENCE_TEMPLATE"], $right_pos + $right_len);

	$vis[$num] = vsprintf("Primer {$num}: %s>%s>%s_%s_%s<%s<%s", $vis[$num]);
}

//blast against transcript database
// Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
$blast_hits = [];
for ($i = 0; $i < $n; ++ $i)
{
	$blast_hits["$i"] = [ "left" => [] , "right" => [] ];
}

list($stdout, ) = $parser->exec("blastn", "-task blastn-short -db /mnt/storage3/users/ahmattj1/GRCh37_cDNA_blast -query {$primers_fasta} -outfmt 6 -evalue 0.1", true);
foreach ($stdout as $line)
{
	$fields = array_combine([ "query", "subject", "identity", "aln_length", "mismatches", "gaps", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"],
		explode("\t", trim($line)));

		//parse read name to determine primer num and orientation
	if (!preg_match("/primer([0-9])+_(left|right)_[0-9]+bp/", $fields["query"], $matches))
	{
		continue;
	}
	$num = $matches[1];
	$which = $matches[2];
	array_push($blast_hits[$num][$which], $fields["subject"]);
	$blast_hits[$num][$which] = array_unique($blast_hits[$num][$which]);
}

//spliced alignment against genomic reference (STAR)
$star_out = $parser->tempFile(".sam");
//$star_out = "test.sam";
$STAR_tmp_folder = $parser->tempFolder();
$genome_dir = "/mnt/storage2/megSAP/data/genomes/STAR/{$build}";
$parser->exec(get_path("STAR"), "--genomeLoad LoadAndKeep --genomeDir {$genome_dir} --outFileNamePrefix {$STAR_tmp_folder}/ --readFilesIn {$primers_fasta} --outSAMtype SAM --outStd SAM > {$star_out}", true);
//$parser->log("Primers SAM:", file($bam));

//create output BAM if requested
if (isset($bam))
{
	$parser->sortBam($star_out, $bam, 1);
	$parser->indexBam($bam, 1);
}

//prepare output table
$out_table = new Matrix();
$out_table->setHeaders([
	"pair",
	"product size",
	"primer",
	"sequence",
	"Tm",
	"GC%",
	"genomic_coordinates",
	"junction",
	"warnings",
	"other cDNA hits (blastn)",
	"visual sequence"
	]);
$out_table->addComment("Primer design: $prefix");
$out_table->addComment("Target transcript: $transcript_id");
$out_table->addComment("Variant: c.{$variant_position}, {$variant_length}bp (cDNA/CDS offset: {$offset})");

//parse alignment
foreach (file($star_out) as $line)
{
	$line = trim($line);

	//skip comments
	if (starts_with($line, "@"))
	{
		continue;
	}

	//BAM fields
	$fields = explode("\t", $line);
	list($readid, $flags, $chr, $pos, $qual, $cigar) = array_slice($fields, 0, 6);
	$tags = array_slice($fields, 11);

	//parse read name to determine primer num and orientation
	if (!preg_match("/primer([0-9])+_(left|right)_[0-9]+bp/", $readid, $matches))
	{
		continue;
	}
	$num = $matches[1];
	$which = $matches[2];

	//calculate alignment positions from CIGAR
	$regions = [];
	$spliced = preg_match_all("/([0-9]+)(M|N)/", $cigar, $matches, PREG_SET_ORDER) > 1;
	$start = $pos;
	foreach ($matches as $component)
	{
		$len = $component[1];
		$end = $start + $len - 1;

		if ($component[2] == "M")
		{
			$regions[] = "{$chr}:{$start}-{$end}";
		}
		
		$start = $end + 1;
	}

	//checks
	$warnings = [];

	//SNPs
	$db_file = get_path("GRCh37_data_folder") . "/dbs/1000G/1000g_v5b.vcf.gz";
	$region_snp = str_replace("chr", "", implode(" ", $regions));
	list($snps, ) = $parser->exec("tabix", "$db_file $region_snp", false);
	if (count($snps) > 0)
	{
		$warnings[] = "SNP";
	}
	
	//number of alignments, mismatches
	foreach ($tags as $tag)
	{
		$tag_fields = explode(":", $tag);
		if ($tag_fields[0] === "NH" && $tag_fields[2] > 1)
		{
			$warnings[] = "multiple alignments";
		}
		if ($tag_fields[0] === "nM" && $tag_fields[2] > 0)
		{
			$warnings[] = "mismatch";
		}
	}
	
	//check strand of forward, reverse primers
	if ($plus_strand && (($which == "left" && $flags & 0x10) || ($which == "right" && !($flags & 0x10)))) {
		$warnings[] = "wrong strand";
	}
	if (!$plus_strand && (($which == "right" && $flags & 0x10) || ($which == "left" && !($flags & 0x10)))) {
		$warnings[] = "wrong strand";
	}

	//check if transcript is in blast hits
	if (!in_array($transcript_id, $blast_hits[$num][$which]))
	{
		$warnings[] = "no blast hit";
	}

	$out_hits = $blast_hits[$num][$which];
	$key = array_search($transcript_id, $out_hits);
	if ($key !== FALSE)
	{
		unset($out_hits[$key]);
	}

	//add to output
	$out_table->addRow([
		$num,
		$results[$num]["PAIR"]["PRODUCT_SIZE"],
		$which,
		$results[$num][strtoupper($which)]["SEQUENCE"],
		$results[$num][strtoupper($which)]["TM"],
		$results[$num][strtoupper($which)]["GC_PERCENT"],
		implode(",", $regions),
		$spliced ? "yes" : "no",
		implode(", ", $warnings),
		implode(",", $out_hits),
		$vis[$num]
		]);
}

$out_table->toTSV($out);

?>