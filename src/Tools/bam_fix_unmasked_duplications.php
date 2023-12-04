<?php
/** 
	@page bam_fix_unmasked_duplications
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("bam_fix_unmasked_duplications", "Fixes BAM files mapped to GRCh38 with unmaked false duplications.");
$parser->addInfile("in", "Input BAM/CRAM file.", false);
$parser->addInfile("in_ref", "Input reference genome used to create the input BAM file.", false);
$parser->addInfile("reg", "Regions of the input BAM which the reads are re-mapped.", true);
$parser->addOutfile("out", "Output CRAM file.", false);
$parser->addInt("threads", "Number of threads to use.", true, 1);
$parser->addFlag("debug", "Debug more. Temp files in output folder and more log output.");
extract($parser->parse($argv));

//init
$ngsbits = get_path("ngs-bits");
$samtools = get_path("samtools");
$bwa = get_path("bwa-mem2");
$out_folder = realpath(dirname($out));
$ref = genome_fasta("GRCh38", true, false);

//set default region
if ($reg=="")
{
	$reg = repository_basedir()."/data/misc/regions_affected_by_GRCh38_false_duplications.bed";
	print "Regions parameter unset - using {$reg}\n";
}

//extract read names
$time_start = microtime(true);
print "extracting read names...\n";
$read_ids = $debug ? $out_folder."/read_ids.txt" : $parser->tempFile("_read_names.txt");
exec2("{$samtools} view -L {$reg} -@ {$threads} -M -T {$in_ref} {$in} | cut -f1 | sort | uniq > {$read_ids}");
list($stdout) = exec2("wc -l < $read_ids");
print "  extracted ".trim(implode("", $stdout))." read names\n";
print "  took ".time_readable(microtime(true)-$time_start)."\n";

//split input BAM according to read names
$time_start = microtime(true);
print "splitting input BAM...\n";
$bam_reg = $debug ? $out_folder."/reg.bam" : $parser->tempFile("_reg.bam");
$bam_other = $debug ? $out_folder."/other.bam" : $parser->tempFile("_other.bam");
$parser->exec("{$ngsbits}/BamExtract", "-in {$in} -ref {$in_ref} -ids {$read_ids} -out {$bam_reg} -out2 {$bam_other}");
print "  took ".time_readable(microtime(true)-$time_start)."\n";

//convert BAM to FASTQ
$time_start = microtime(true);
print "converting BAM to FASTQ\n";
$fastq1 = $debug ? $out_folder."/R1.fastq.gz" : $parser->tempFile("_R1.fastq.gz");
$fastq2 = $debug ? $out_folder."/R2.fastq.gz" : $parser->tempFile("_R2.fastq.gz");
$parser->exec("{$ngsbits}/BamToFastq", "-in {$bam_reg} -out1 {$fastq1} -out2 {$fastq2}");
print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
//re-map extracted reads
$time_start = microtime(true);
print "re-mapping reads\n";
$bam_mapped = $debug ? $out_folder."/mapped.bam" : $parser->tempFile("_mapped.bam");
$args = [];
$args[] = "-in1 {$fastq1}";
$args[] = "-in2 {$fastq2}";
$args[] = "-out {$bam_mapped}";
$args[] = "-sample ".basename2($in);
$args[] = "-threads {$threads}";
$args[] = "-dedup";
$parser->execTool("NGS/mapping_bwa.php", implode(" ", $args));
print "  took ".time_readable(microtime(true)-$time_start)."\n";

//merge sorted BAMs
$time_start = microtime(true);
print "merging BAMs\n";
exec2("{$samtools} merge --write-index -@ {$threads} -f -h {$bam_mapped} --reference {$ref} -o {$out} {$bam_mapped} {$bam_other}");
print "  took ".time_readable(microtime(true)-$time_start)."\n";

?>

