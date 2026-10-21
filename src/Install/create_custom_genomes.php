<?php

/**
  @page create_custom_genomes
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("create_custom_genomes", "Creates custom masked genomes for targeted calling.");
$parser->addInfile("out", "Output folder to store genomes.", false, false);
$parser->addInfile("target", "BED file with target region for variant calling (Region to keep).", false);

//optional
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));


//split target and source regions by gene
function split_bed_by_gene($bed_file)
{
	$region_by_gene = [];
	$buffer = file($bed_file, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
	foreach ($buffer as $line) 
	{
		list(, , , $gene) = explode("\t", $line);
		if (!isset($gene)) $region_by_gene[$gene] = [];
		$region_by_gene[$gene][] = $line;
	}

	return $region_by_gene;
}

function get_chr_by_gene($target_by_gene)
{
	$chr_by_gene = [];
	foreach ($target_by_gene as $gene => $roi) 
	{
		if (!isset($chr_by_gene[$gene])) $chr_by_gene[$gene] = [];
		foreach ($roi as $line)
		{
			list($chr, , ) = explode("\t", $line);
			$chr_by_gene[$gene][$chr] = true;
		}
	}
	foreach ($chr_by_gene as $gene => $chrs) 
	{
		
		$tmp = array_keys($chrs);
		sort($tmp);
		$chr_by_gene[$gene] = $tmp;
	}
	return $chr_by_gene;
}

//create masked genome for each gene
function get_mask_genome($chromosomes, $region2mask, $gene)
{
	global $parser;
	global $genome;
	global $build;
	global $out;
	print "generating custom genome for gene {$gene} (".implode(",", $chromosomes).")...\n";
	//slice genome
	$tmp_genome = $parser->tempFile(".fa");
	$parser->execApptainer("samtools", "samtools faidx", "-o {$tmp_genome} {$genome} ".implode(" ", $chromosomes), [$genome]);
	//mask genome slice
	$masked_genome = "{$out}/{$gene}_{$build}.fa";
	$parser->execApptainer("ngs-bits", "FastaMask", "-in {$tmp_genome} -reg {$region2mask} -out {$masked_genome}", [], [dirname($masked_genome)]);
	$parser->execApptainer("samtools", "samtools faidx", $masked_genome, [], [dirname($masked_genome)]);
	//create BWA index
	$parser->execApptainer("bwa-mem2", "bwa-mem2".trim(get_path("bwa_mem2_suffix", false)), " index {$masked_genome}", [], [dirname($masked_genome)]);
}

//init
$genome = genome_fasta($build);
$chr_size = genome_chr_sizes($genome);

// create output folder
if (!file_exists($out))
{
	exec2("mkdir -p $out");
}

//split target/source by gene
$target_by_gene = split_bed_by_gene($target);

//check for missing/mismatching genes
$genes_target = array_keys($target_by_gene);
//get covered chrs per gene
$chr_by_gene = get_chr_by_gene($target_by_gene);

//create custom genomes
$custom_genomes = [];
foreach ($target_by_gene as $gene => $roi) 
{
	//store region to keep
	$tmp_bed = $parser->tempFile(".bed");
	file_put_contents($tmp_bed, implode("\n", $roi));
	$tmp_chr_bed = $parser->tempFile(".bed");
	$chr_region = [];
	foreach ($chr_by_gene[$gene] as $chr) 
	{
		$chr_region[] = "{$chr}\t0\t".$chr_size[$chr];
	}
	file_put_contents($tmp_chr_bed, implode("\n", $chr_region));
	$region2mask = $parser->tempFile(".bed");
	$parser->execApptainer("ngs-bits", "BedSubtract", "-in {$tmp_chr_bed} -in2 {$tmp_bed} -out {$region2mask}", [], []);

	get_mask_genome($chr_by_gene[$gene], $region2mask, $gene);
}

?>
