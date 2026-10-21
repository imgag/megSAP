<?php

/**
  @page vc_mrjd
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_mrjd", "Remap homologous regions and call variants in these regions.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addOutfile("out", "Output VCF.GZ file.", false);

$parser->addInfile("target", "BED file with target region for variant calling.", false);
$parser->addInfile("source", "BED file with source region for remapping calling.", false);
$parser->addInfile("masked_genome_folder", "Folder containing custom masked genomes for each gene (requires bwa-mem2 index).", false);

//optional
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addInt("threads", "Number of threads used.", true, 1);
$parser->addFlag("keep_bam", "Store remapped BAM file in sample folder.");
$parser->addInt("max_mq", "Maximum mapping quality used to filter reads in `-source`", true, 50);
$parser->addInfile("blacklist", "VCF which will be removed from result.", true);
$parser->addFlag("flag_blacklist", "Flag variants from the blacklist instead of removing them.");
extract($parser->parse($argv));


//split target and source regions by gene
function split_bed_by_gene($bed_file)
{
	$region_by_gene = [];
	$buffer = file($bed_file, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
	foreach ($buffer as $line) 
	{
		$parts = explode("\t", $line);
		if (count($parts) < 4) trigger_error("Line '{$line}' does not contain gene info!", E_USER_ERROR);
		$gene = $parts[3];
		if (!isset($region_by_gene[$gene])) $region_by_gene[$gene] = [];
		$region_by_gene[$gene][] = $line;
	}

	return $region_by_gene;
}


function extract_reads($roi_source, $roi_target)
{
	global $parser;
	global $in;
	global $genome;
	global $threads;
	global $max_mq;

	if (count($roi_source) < 1) 
	{
		trigger_error("Source region is empty!", E_USER_WARNING);
		$source_reads_bam_filtered = $parser->tempFile("_source_reads_filtered.bam");
		$parser->execApptainer("samtools", "samtools view", "-T {$genome} -H -o {$source_reads_bam_filtered} {$in}", [$in, $genome]);
	}
	else
	{
		// extract low quality reads (read names) in source region
		$source_reads_bam = $parser->tempFile("_source_reads.bam");
		$tmp_source_bed = $parser->tempFile(".bed");
		file_put_contents($tmp_source_bed, implode("\n", $roi_source));
		$parser->execApptainer("samtools", "samtools view", "-T {$genome} -h --region-file {$tmp_source_bed} -@ {$threads} -o {$source_reads_bam} {$in}", [$in, $genome]);
		$source_reads_bam_filtered = $parser->tempFile("_source_reads_filtered.bam");
		$parser->execApptainer("ngs-bits", "BamFilter", "-in {$source_reads_bam} -out {$source_reads_bam_filtered} -maxMQ {$max_mq}");
	}
	

	//extract all reads in target region
	$target_reads_bam = $parser->tempFile("_target_reads.bam");
	$tmp_target_bed = $parser->tempFile(".bed");
	file_put_contents($tmp_target_bed, implode("\n", $roi_target));
	$parser->execApptainer("samtools", "samtools view", "-T {$genome} -h --region-file {$tmp_target_bed} -@ {$threads} -o {$target_reads_bam} {$in}", [$in, $genome]);

	//merge in one file and convert to fastq
	$combined_bam = $parser->tempFile("_combined.bam");
	$parser->execApptainer("samtools", "samtools cat", "-o {$combined_bam} {$source_reads_bam_filtered} {$target_reads_bam}");
	$reads_r1 = $parser->tempFile("_R1.fastq.gz");
	$reads_r2 = $parser->tempFile("_R2.fastq.gz");
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$combined_bam} -out1 {$reads_r1} -out2 {$reads_r2}");

	return [$reads_r1, $reads_r2];
}

function remap_reads($reads_r1, $reads_r2, $masked_genome)
{
	global $parser;
	global $threads;
	global $build;
	//remap to modified genome
	$pipeline = [];
	$bwa_params = "mem {$masked_genome} -K 100000000 -Y -t {$threads} -v 2";
	$suffix = trim(get_path("bwa_mem2_suffix", false));

	$pipeline[] = ["", $parser->execApptainer("bwa-mem2", "bwa-mem2{$suffix}", "{$bwa_params} {$reads_r1} {$reads_r2}", [$masked_genome], [], true)];
	$pipeline[] = ["", $parser->execApptainer("samblaster", "samblaster", "", [], [], true)]; 
	$tmp_unsorted = $parser->tempFile("_unsorted.bam");
	$pipeline[] = ["", $parser->execApptainer("samtools", "samtools", "view -Sb - > $tmp_unsorted", [], [], true)];
	$parser->execPipeline($pipeline, "remapping");

	$remapped_bam = $parser->tempFile(".bam");
	$parser->sortBam($tmp_unsorted, $remapped_bam, $threads, $build);
	$parser->indexBam($remapped_bam, $threads);

	return $remapped_bam;
}

function vc_calling($bam, $roi_target)
{
	global $parser;
	global $threads;
	global $build;

	//perform variant calling
	$tmp_bed = $parser->tempFile("_target_region.bed");
	file_put_contents($tmp_bed, implode("\n", $roi_target));
	$tmp_vcf = $parser->tempFile("_targeted.vcf.gz");
						
	$args = [];
	$args[] = "-bam {$bam}";
	$args[] = "-out {$tmp_vcf}";
	$args[] = "-build {$build}";
	$args[] = "-threads {$threads}";
	$args[] = "-target {$tmp_bed}";
	$args[] = "-min_af 0.05";
	$args[] = "-min_mq 0";
	$args[] = "-min_qsum 0";
	$args[] = "-min_bq 10";
	$args[] = "-no_ploidy";
	$parser->execTool("Tools/vc_freebayes.php", implode(" ", $args));

	return $tmp_vcf;
}

function remove_blacklist_variants($vcf, $blacklist, $flag_blacklist)
{
	global $parser;
	global $threads;

	$filtered_vcf = $parser->tempFile("_targeted.vcf.gz");
		
	if ($flag_blacklist)
	{
		//only flag
		if (ends_with($blacklist, ".vcf.gz")) $blacklist_vcf = $blacklist;
		else
		{
			$blacklist_vcf = $parser->tempFile("_blacklist.vcf.gz");
			$parser->exec("bgzip", "-@ {$threads} -l 7 -c {$blacklist} > {$blacklist_vcf}");
			$parser->execApptainer("htslib", "tabix", "-f -p vcf {$blacklist_vcf }");
		}
	
		$pipeline = [];
		$pipeline[] = ["zcat", $vcf];
		$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfAnnotateFromVcf", 
		"-source {$blacklist_vcf} -existence_only -existence_key_name HIGH_POP_AF -threads {$threads}", [$blacklist_vcf], [], true)];
		$pipeline[] = ["bgzip", "-l 7 -c > {$filtered_vcf}"];
		$parser->execPipeline($pipeline, "blacklist flagging");

	}
	else
	{
		//remove
		if (ends_with($blacklist, ".vcf.gz"))
		{
			$blacklist_vcf = $parser->tempFile("_blacklist.vcf");
			$parser->exec("zcat", "{$blacklist} > {$blacklist_vcf}");
		}
		else $blacklist_vcf = $blacklist;
		$pipeline = [];
		$pipeline[] = ["zcat", $vcf];
		$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfSubtract", "-in2 {$blacklist_vcf}", [$blacklist_vcf], [], true)];
		$pipeline[] = ["bgzip", "-l 7 -c > {$filtered_vcf}"];
		$parser->execPipeline($pipeline, "blacklist filtering");
	}
	
	return $filtered_vcf;
}

//init
$genome = genome_fasta($build);

//split target/source by gene
$target_by_gene = split_bed_by_gene($target);
$source_by_gene = split_bed_by_gene($source);
//add missing source keys
foreach ($target_by_gene as $gene => $roi) if (!isset($source_by_gene[$gene])) $source_by_gene[$gene] = [];

//check for missing/mismatching genes
$genes_target = array_keys($target_by_gene);
$genes_source = array_keys($source_by_gene);
sort($genes_target);
sort($genes_source);
if($genes_target != $genes_source) 
{
	trigger_error("Gene list of target and source region mismatch:\n\ttarget:\t".implode(",", $genes_target)."\n\tsource:\t".implode(",", $genes_source), E_USER_ERROR);
}

//extract reads
$extracted_reads = [];
foreach ($source_by_gene as $gene => $roi_source) 
{
	print "extracting reads for gene {$gene}...\n";
	$roi_target = $target_by_gene[$gene];
	$extracted_reads[$gene] = extract_reads($roi_source, $roi_target);
}

//remap reads
$remapped_bams = [];
foreach ($source_by_gene as $gene => $roi_source) 
{
	print "re-mapping reads for gene {$gene}...\n";
	$reads_r1 = $extracted_reads[$gene][0];
	$reads_r2 = $extracted_reads[$gene][1];
	$custom_genome = "{$masked_genome_folder}/{$gene}_{$build}.fa";
	if (!file_exists($custom_genome)) trigger_error("Custom genome file '{$custom_genome}' does not exist!", E_USER_ERROR);
	$remapped_bams[$gene] = remap_reads($reads_r1, $reads_r2, $custom_genome);
}

//vc calling
$vcfs = [];
foreach ($target_by_gene as $gene => $roi_target) 
{
	print "recalling gene {$gene}...\n";
	$vcfs[$gene] = vc_calling($remapped_bams[$gene], $roi_target);
}


//merge VCF
$tmp_vcf = $parser->tempFile(".vcf.gz");
$parser->execApptainer("ngs-bits", "VcfAdd", "-skip_duplicates -out {$tmp_vcf} -in ".implode(" ", $vcfs), [], []);
$parser->execApptainer("ngs-bits", "VcfSort", "-compression_level 7 -in {$tmp_vcf} -out {$tmp_vcf}", [], []);

//remove blacklist variants
if (!empty($blacklist))
{
	$tmp_vcf = remove_blacklist_variants($tmp_vcf, $blacklist, $flag_blacklist);
}

$parser->moveFile($tmp_vcf, $out);
$parser->execApptainer("htslib", "tabix", "-f -p vcf {$out}", [], [dirname($out)]);

if ($keep_bam)
{
	$out_cram = dirname($out)."/".basename($out, ".vcf.gz").".cram";
	$parser->execApptainer("samtools", "samtools merge", "-f --reference {$genome} -@ {$threads} -O cram,version=3.0 -o {$out_cram} ".implode(" ", $remapped_bams), [$genome], [dirname($out_cram)]);
	$parser->execApptainer("samtools", "samtools index", "-@ {$threads} {$out_cram}", [], [dirname($out_cram)]);
}


?>
