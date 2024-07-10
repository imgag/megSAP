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
$parser->addFlag("remap_missing_chr", "Also remap all reads on chromosomes which are not present in the target reference.");
$parser->addFlag("replace_readgroup", "Replace the read group. Either using '-sample' or basename of input BAM file.");
$parser->addString("sample", "Sample name used as read group. If empty name of input bam is used", true, "");
$parser->addFlag("debug", "Debug more. Temp files in output folder and more log output.");
extract($parser->parse($argv));

//init
$ngsbits = get_path("ngs-bits");
$samtools = get_path("samtools");
$bwa = get_path("bwa-mem2");
$out_folder = realpath(dirname($out));
$ref = genome_fasta("GRCh38", true, false);
if ($sample == "") $sample = basename2($in);
print $sample;
$removed_chromosomes = array();

//set default region
if ($reg=="")
{
	$reg = repository_basedir()."/data/misc/regions_affected_by_GRCh38_false_duplications.bed";
	print "Regions parameter unset - using {$reg}\n";
}

if ($remap_missing_chr)
{
	//get missing chromosomes
	$tmp = $parser->tempFile(".tsv");
	list($stdout, $stderr, $exit_code) = exec2("diff {$in_ref}.fai {$ref}.fai | egrep '^<' | cut -d ' ' -f2 | cut -f1-2 > {$tmp}", false);
	if ($exit_code > 1 || $exit_code < 0)
	{
		trigger_error("During extraction of missing chromosomes the following error occured: \n".implode("\n", $stderr), E_USER_ERROR);
	}
	elseif ($exit_code == 1)
	{
		//generate BED file
		$missing_chr_content = file($tmp);
		$bed_content = array();
		foreach ($missing_chr_content as $line) 
		{
			$split_line = explode("\t", $line);
			$chr = trim($split_line[0]);
			$start = 0;
			$end = intval($split_line[1]);
			$bed_content[] = implode("\t", array($chr, $start, $end));
			$removed_chromosomes[] = $chr;
		}
		$tmp_bed_missing_chr = $parser->tempFile(".bed");
		file_put_contents($tmp_bed_missing_chr, implode("\n", $bed_content));
		
		//merge BED files
		$merged_bed = $parser->tempFile(".bed");
		$parser->exec(get_path("ngs-bits")."BedAdd", "-in {$reg} {$tmp_bed_missing_chr} -out {$merged_bed}");

		//set merged bed as input
		$reg = $merged_bed;
	}
	else
	{
		trigger_error("No missing chromosomes found.", E_USER_NOTICE);
	}

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
$args[] = "-sample {$sample}";
$args[] = "-threads {$threads}";
$args[] = "-dedup";
$parser->execTool("NGS/mapping_bwa.php", implode(" ", $args));
print "  took ".time_readable(microtime(true)-$time_start)."\n";

//determine remapped regions:
$time_start = microtime(true);
$remapped_region = $out_folder."/remapped_regions.bed";
print "get region of remapped reads\n";
$parser->exec("{$ngsbits}/BedHighCoverage", "-bam {$bam_mapped} -out {$remapped_region} -cutoff 1 -threads {$threads}");
print "  took ".time_readable(microtime(true)-$time_start)."\n";

$merged_bam = $debug ? $out_folder."/merged.bam" : $parser->tempFile("_mapped.bam");
if (!$replace_readgroup) $merged_bam = $out; //write directly to output folder

//merge sorted BAMs
$time_start = microtime(true);
print "merging BAMs\n";
exec2("{$samtools} merge --write-index -@ {$threads} -f -h {$bam_mapped} --reference {$ref} -o {$merged_bam} {$bam_mapped} {$bam_other}");
print "  took ".time_readable(microtime(true)-$time_start)."\n";

if ($replace_readgroup)
{
	print "replace RG\n";
	list($in_header) = exec2("{$samtools} view -H {$merged_bam}");
	
	//extract @SQ lines from mapped file
	list($in_header2) = exec2("{$samtools} view -H {$bam_mapped}");
	$new_sq_lines = array();
	foreach ($in_header2 as $line) 
	{
		if (starts_with($line, "@SQ")) $new_sq_lines[] = $line;
	}
	//remove unused header lines
	$out_header = array();
	$sq_replaced = false;
	foreach ($in_header as $line) 
	{
		if (!starts_with($line, "@RG") && !starts_with($line, "@SQ")) $out_header[] = $line; //keep all non read group and chr lines
		//keep only remaining chr
		if (!$sq_replaced && starts_with($line, "@SQ"))
		{
			$out_header = array_merge($out_header, $new_sq_lines);
			$sq_replaced = true;
		}
		//keep only read group entry of given sample
		if (starts_with($line, "@RG\tID:{$sample}\t")) $out_header[] = $line;		 
	} 
	$tmp_header = $parser->tempFile("_header.sam");
	file_put_contents($tmp_header, implode("\n", $out_header));
	
	//re-header & replace RG of all reads
	$tmp_fixed_bam = $parser->tempFile("_fixed_rg.bam");
	$parser->exec($samtools, "addreplacerg -@ {$threads} -R {$sample} -m overwrite_all -o {$tmp_fixed_bam} {$merged_bam}");
	$tmp_fixed_bam2 = $debug ? $out_folder."/fixed_rg_header.bam" : $parser->tempFile("_fixed_header.bam");
	$parser->exec($samtools, "reheader {$tmp_header} {$tmp_fixed_bam} > {$tmp_fixed_bam2}");
	$parser->exec($samtools, "view -h --write-index -@ {$threads} -T {$ref} -C -o {$out} {$tmp_fixed_bam2}");
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
}

?>

