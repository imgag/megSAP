<?php

/**
 * vc_viral_load
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_viral_load", "Estimate viral load from aligned data by aligning to viral sequence reference.");
$parser->addInfile("in",  "Input file in BAM format.", false);
$parser->addOutfile("viral_bam", "Output file in BAM format.", false);
$parser->addOutFile("viral_bam_raw", "Undeduplicated output in BAM format in case of barcode correction.", true);
$parser->addOutfile("viral_cov", "Output file in TSV format.", false);

$parser->addInfile("in_qcml", "Mapping statistics of input data in qcML format for relative coverage values.", true);
$parser->addFlag("barcode_correction", "Run UMI-specific deduplication steps.");
$parser->addStringArray("viral_chrs", "Viral chromosome names or region specifiers in original alignment, i.e. for HHV-4.", true, ["chrNC_007605"]);

$parser->addString("build_viral", "Build name of viral references.", true, "somatic_viral");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));

if ($barcode_correction && !isset($viral_bam_raw))
{
    $viral_bam_raw = $parser->tempFile("_viral_before_dedup.bam");
}

//data setup
$parser->execTool("Tools/data_setup.php", "-build {$build_viral}");

//target file for viral sequences
$viral_enrichment = get_path("data_folder") . "/enrichment/{$build_viral}.bed";

//extract unaligned reads
$filtered_bam1 = $parser->tempFile("_filtered1.bam");
$parser->exec(get_path("samtools"), "view -u -U {$filtered_bam1} -F 12 {$in} > /dev/null", true);

//extract alignments from viral chromosomes/regions, i.e. chrNC_007605
if (count($viral_chrs) > 0)
{
    $filtered_bam2 = $parser->tempFile("_filtered2.bam");
    $viral_regions = implode(" ", $viral_chrs);
    $parser->exec(get_path("samtools"), "view -u -o {$filtered_bam2} {$in} {$viral_regions}", true);
    
    $filtered_bam = $parser->tempFile("_filtered.bam");
    $parser->exec(get_path("samtools"), "merge -u {$filtered_bam} {$filtered_bam1} {$filtered_bam2}", true);
}
else {
    $filtered_bam = $filtered_bam1;
}

//generate FASTQ
$filtered_r1 = $parser->tempFile("_filtered_R1.fastq.gz");
$filtered_r2 = $parser->tempFile("_filtered_R2.fastq.gz");
list($stdout, ) =$parser->exec(get_path("ngs-bits")."/BamToFastq", "-in {$filtered_bam} -out1 {$filtered_r1} -out2 {$filtered_r2}", true);

//align to viral sequences
$viral_tmp1 = $parser->tempFile("_viral_tmp.bam");
$mapping_dedup = $barcode_correction ? "" : "-dedup";
$parser->execTool("NGS/mapping_bwa.php", "-in1 {$filtered_r1} -in2 {$filtered_r2} -out {$viral_tmp1} {$mapping_dedup} -build {$build_viral} -threads {$threads}");

//run UMI-based deduplication if requested
if ($barcode_correction)
{
    $viral_tmp2 = $parser->tempFile("_viral_tmp.bam");
    $parser->exec(get_path("python3")." ".get_path("umiVar2")."/barcode_correction.py", "--infile $viral_tmp1 --outfile $viral_tmp2", true);
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
$parser->exec(get_path("ngs-bits")."BamClipOverlap", "-in $viral_tmp2 -out $viral_tmp3 -overlap_mismatch_basen", true);
$parser->sortBam($viral_tmp3, $viral_bam, $threads);
$parser->indexBam($viral_bam, $threads);

//run variant calling
$viral_vcfgz = $parser->tempFile("_viral.vcf.gz");
$parser->execTool("NGS/vc_freebayes.php", "-bam {$viral_bam} -target {$viral_enrichment} -build {$build_viral} -no_ploidy -min_af 0.01 -out {$viral_vcfgz}");
$viral_vcf = $parser->tempFile("_viral.vcf");
$parser->exec("gzip", "-cd {$viral_vcfgz} > {$viral_vcf}", true);

//create output report file
$viral_cov_tmp = $parser->tempFile("_viral_cov.bed");
$parser->exec(get_path("ngs-bits")."/BedCoverage", "-out {$viral_cov_tmp} -min_mapq 30 -decimals 4 -in {$viral_enrichment} -bam {$viral_bam} -threads {$threads}", true);

$viral_reads_tmp = $parser->tempFile("_viral_reads.bed");
$parser->exec(get_path("ngs-bits")."/BedReadCount", "-out {$viral_reads_tmp} -min_mapq 30 -in {$viral_enrichment} -bam {$viral_bam}", true);

$viral_result_tmp = $parser->tempFile("_viral_result.bed");
$pipeline = [];
$pipeline[] = [ get_path("ngs-bits")."/BedAnnotateFromBed", "-in {$viral_enrichment} -no_duplicates -in2 {$viral_reads_tmp}" ];
$pipeline[] = [ get_path("ngs-bits")."/BedAnnotateFromBed", "-in2 {$viral_cov_tmp} -no_duplicates -col 5 -out {$viral_result_tmp}" ];
$parser->execPipeline($pipeline, "result file");

//extract average target coverage from provided qcML file
$avg_target_cov = 0.0;
if (isset($in_qcml)) {
    $sxml = simplexml_load_file($in_qcml);
    foreach ($sxml->runQuality->qualityParameter as $qp)
    {
        if ($qp["accession"] == "QC:2000025")
        {
            $avg_target_cov = floatval($qp["value"]);
        }
    }
}

//add additional columns with relative coverage and name
$handle_out = fopen2($viral_cov, "w");
preg_match('/.+ : ([0-9]+)/', $stdout[0], $matches);
fwrite($handle_out,         "##number of reads used for analysis:        " . $matches[1] . "\n");
fwrite($handle_out, sprintf("##average target region depth of input BAM: %.4f\n", $avg_target_cov));
fwrite($handle_out, "#chr\tstart\tend\tname\treads\tcoverage\tcoverage_rel\tmismatches\tidentity%\n");

$handle = fopen2($viral_result_tmp, "r");
while ($line = fgets($handle))
{
    if (starts_with($line, "#chr")) continue;
    $fields = explode("\t", nl_trim($line));

    // relative coverage
    $fields[] = sprintf("%.4f", $fields[5] / $avg_target_cov);

    // mismatches
    $region = $fields[0] . ":" . ($fields[1]+1) . "-" .  $fields[2];
    list($stdout) = $parser->exec(get_path("ngs-bits")."VcfFilter", "-in {$viral_vcf} -reg '{$region}'", false);
    $variant_count = 0;
	foreach($stdout as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		++$variant_count;
	}
    $fields[] = $variant_count;
    $fields[] = 100 * (1 - $variant_count / ($fields[2] - $fields[1]));
    fwrite($handle_out, implode("\t", $fields) . "\n");
}
fclose($handle);
fclose($handle_out);

?>