<?php
/** 
	@page index_genome
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("index_genome", "Indexes a genome FASTA file with BWA and samtools.");
$parser->addInfile("in", "FASTA file to index.", false);
$parser->addFlag("mask", "Mask false duplications in GRCh38 (see https://www.nature.com/articles/s41587-021-01158-1).");
extract($parser->parse($argv));

if ($mask)
{
	$exclusion_bed = $parser->tempFile("_exclusion.bed");
	exec2("wget -O {$exclusion_bed} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_GRC_exclusions.bed");
	$tmp = $parser->tempFile("_masked.fa");
	exec2(get_path("bedtools")." maskfasta -fi {$in} -bed {$exclusion_bed} -fullHeader -fo {$tmp}");
	$parser->moveFile($tmp, $in);
}

if (get_path("use_bwa1"))
{
	exec2(get_path("bwa")." index -a bwtsw {$in}");
}
else
{
	exec2(get_path("bwa-mem2")." index {$in}");
}

exec2(get_path("samtools")." faidx {$in}");
exec2("md5sum -b {$in} > {$in}.md5");

?>
