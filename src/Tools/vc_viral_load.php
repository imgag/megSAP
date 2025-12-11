<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_viral_load", "Estimate viral load from aligned data by aligning to viral sequence reference.");
$parser->addInfile("in",  "Input file in BAM format.", false);
$parser->addOutfile("viral_bam", "Output file in BAM format.", false);
$parser->addOutFile("viral_bam_raw", "Undeduplicated output in BAM format in case of barcode correction.", true);
$parser->addOutfile("viral_cov", "Output file in TSV format.", false);
$parser->addFloat("avg_target_cov", "Average read depth of tumor BAM file.", false);
//optional
$parser->addFlag("barcode_correction", "Run UMI-specific deduplication steps.");
$parser->addStringArray("viral_chrs", "Viral chromosome names or region specifiers in original alignment, i.e. for HHV-4.", true, ["chrNC_007605"]);
$parser->addString("build_viral", "Build name of viral references.", true, "somatic_viral");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));


if ($barcode_correction && !isset($viral_bam_raw))
{
    $viral_bam_raw = $parser->tempFile("_viral_before_dedup.bam");
}

//create viral genome index if missing
$viral_genome  = get_path("data_folder")."/genomes/{$build_viral}.fa";
if(!file_exists($viral_genome)) trigger_error("Viral reference genome mnissing: {$viral_genome}.", E_USER_ERROR);
if (file_exists($viral_genome) && !file_exists($viral_genome.".bwt.2bit.64"))
{
	$parser->execTool("Install/index_genome.php", "-in {$viral_genome}");
}

//target file for viral sequences
$viral_enrichment = get_path("data_folder") . "/enrichment/{$build_viral}.bed";
if(!file_exists($viral_enrichment)) trigger_error("Viral target region missing: {$viral_enrichment}.", E_USER_ERROR);

//extract unaligned reads
$genome = genome_fasta("GRCh38");
$filtered_bam1 = $parser->tempFile("_filtered1.bam");
$parser->execApptainer("samtools", "samtools view", "-T {$genome} -u -U {$filtered_bam1} -F 12 {$in} > /dev/null", [$genome, $in]);

//extract alignments from viral chromosomes/regions, i.e. chrNC_007605
if (count($viral_chrs) > 0)
{
    $filtered_bam2 = $parser->tempFile("_filtered2.bam");
    $viral_regions = implode(" ", $viral_chrs);
    $parser->execApptainer("samtools", "samtools view", "-T {$genome} -u -o {$filtered_bam2} {$in} {$viral_regions}", [$genome, $in]);
    
    $filtered_bam = $parser->tempFile("_filtered.bam");
    $parser->execApptainer("samtools", "samtools merge", "-u {$filtered_bam} {$filtered_bam1} {$filtered_bam2}");
}
else
{
    $filtered_bam = $filtered_bam1;
}

//generate FASTQ
$filtered_r1 = $parser->tempFile("_filtered_R1.fastq.gz");
$filtered_r2 = $parser->tempFile("_filtered_R2.fastq.gz");
list($stdout, ) = $parser->execApptainer("ngs-bits", "BamToFastq", "-in {$filtered_bam} -out1 {$filtered_r1} -out2 {$filtered_r2}");

//align to viral sequences
$viral_tmp1 = $parser->tempFile("_viral_tmp.bam");
$mapping_dedup = $barcode_correction ? "" : "-dedup";
$parser->execTool("Tools/mapping_bwa.php", "-in1 {$filtered_r1} -in2 {$filtered_r2} -out {$viral_tmp1} {$mapping_dedup} -build {$build_viral} -threads {$threads}");

//run UMI-based deduplication if requested
if ($barcode_correction)
{
    $viral_tmp2 = $parser->tempFile("_viral_tmp.bam");
    $parser->execApptainer("umiVar", "barcode_correction.py", "--infile $viral_tmp1 --outfile $viral_tmp2");
    $parser->indexBam($viral_tmp2, $threads);
    
    $parser->moveFile($viral_tmp1, $viral_bam_raw);
    $parser->moveFile($viral_tmp1 . ".bai", $viral_bam_raw . ".bai");
}
else
{
    $viral_tmp2 = $viral_tmp1;
}

//run BamClipOverlap
$viral_tmp3 = $parser->tempFile("_viral_clipoverlap_tmp.bam");
$parser->execApptainer("ngs-bits", "BamClipOverlap", "-in $viral_tmp2 -out $viral_tmp3 -overlap_mismatch_basen");
$parser->sortBam($viral_tmp3, $viral_bam, $threads, $build_viral);
$parser->indexBam($viral_bam, $threads);

//run variant calling
$viral_vcfgz = $parser->tempFile("_viral.vcf.gz");
$parser->execTool("Tools/vc_freebayes.php", "-bam {$viral_bam} -target {$viral_enrichment} -build {$build_viral} -no_ploidy -min_af 0.01 -out {$viral_vcfgz}");
$viral_vcf = $parser->tempFile("_viral.vcf");
$parser->exec("gzip", "-cd {$viral_vcfgz} > {$viral_vcf}", true);

//create output report file
$viral_cov_tmp = $parser->tempFile("_viral_cov.bed");
$parser->execApptainer("ngs-bits", "BedCoverage", "-out {$viral_cov_tmp} -min_mapq 30 -decimals 4 -in {$viral_enrichment} -bam {$viral_bam} -threads {$threads}", [$viral_enrichment, $viral_bam]);

$viral_reads_tmp = $parser->tempFile("_viral_reads.bed");
$parser->execApptainer("ngs-bits", "BedReadCount", "-out {$viral_reads_tmp} -min_mapq 30 -in {$viral_enrichment} -bam {$viral_bam}", [$viral_enrichment, $viral_bam]);

$viral_result_tmp = $parser->tempFile("_viral_result.bed");
$pipeline = [];
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$viral_enrichment} -no_duplicates -in2 {$viral_reads_tmp}", [$viral_enrichment], [], true)];
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in2 {$viral_cov_tmp} -no_duplicates -col 5 -out {$viral_result_tmp}", [], [], true)];
$parser->execPipeline($pipeline, "result file");

//add additional columns with relative coverage and name
$handle_out = fopen2($viral_cov, "w");
preg_match('/.+ : ([0-9]+)/', $stdout[0], $matches);
fwrite($handle_out, "##number of reads used for analysis: " . $matches[1] . "\n");
fwrite($handle_out, "##average target region depth of input BAM: ".number_format($avg_target_cov, 4)."\n");
fwrite($handle_out, "#chr\tstart\tend\tname\treads\tcoverage\tcoverage_rel\tmismatches\tidentity%\n");

$handle = fopen2($viral_result_tmp, "r");
while ($line = fgets($handle))
{
    if (starts_with($line, "#chr")) continue;
    $fields = explode("\t", nl_trim($line));

    // relative coverage
    $fields[] = number_format($fields[5]/$avg_target_cov, 4);

    // mismatches
    $region = $fields[0] . ":" . ($fields[1]+1) . "-" .  $fields[2];
    list($stdout) = $parser->execApptainer("ngs-bits", "VcfFilter", "-in {$viral_vcf} -reg '{$region}' -ref $genome", [$genome], [], false, false);
    $variant_count = 0;
	foreach($stdout as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		++$variant_count;
	}
    $fields[] = $variant_count;
    $fields[] = number_format(100 * (1 - $variant_count / ($fields[2] - $fields[1])),2);
    fwrite($handle_out, implode("\t", $fields) . "\n");
}
fclose($handle);
fclose($handle_out);

?>