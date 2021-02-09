<?php

/**
 * @page mapping_bismark
 */
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("mapping_bismark", "Bisulfite treated DNA read mapping with Bismark.");

//mandatory arguments
$parser->addString("folder", "Sample folder with FASTQ files.", false);

//optional arguments
$parser->addString("out_folder", "Folder where results are stored, defaults to sample folder.", true, "");
$parser->addString("name", "Sample name, extracted from folder if not specified.", true, "");
$parser->addInfile("system", "Processing system INI file (determined from folder name by default).", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);

extract($parser->parse($argv));

//resolve output folder
if ($out_folder === "")
{
	$out_folder = $folder;
}

//input FASTQ files
$in_for = glob($folder."/*_R1_001.fastq.gz");
$in_rev = glob($folder."/*_R2_001.fastq.gz");

if ((count($in_for) == 0) && (count($in_rev) == 0)) {
	trigger_error("No FASTQ files found!", E_USER_ERROR);
}
if ((count($in_for) != count($in_rev)) && count($in_rev)>0) {
	trigger_error("Mismatching number of R1 and R2 files!", E_USER_ERROR);
}

//use sample id, if possible
$p = basename(realpath($out_folder));
if ($name === "" && preg_match('/^Sample_(.+)/', $p, $matches)) {
	$name = $matches[1];
}
else if ($name === "")
{
	trigger_error("Could not extract sample ID! Use -name to specify." , E_USER_ERROR);
}

$sys = load_system($system, $name);

$in_for_s = implode(" ", $in_for);
$in_rev_s = implode(" ", $in_rev);
$paired = count($in_rev) > 0;

$basename = $out_folder."/".$name;
$qc_fastq = $basename."_stats_fastq.qcML";

//check that adapters are specified
if ($sys["adapter1_p5"]=="" || $sys["adapter2_p7"]=="")
{
	trigger_error("No forward and/or reverse adapter sequence given!\nForward: ".$sys["adapter1_p5"]."\nReverse: ".$sys["adapter2_p7"], E_USER_ERROR);
}

//adapter trimming + QC (SeqPurge for paired-end, ReadQC+skewer for single-end)
if ($paired)
{
	$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
	$fastq_trimmed2 = $parser->tempFile("_trimmed.fastq.gz");
	$seqpurge_params = array(
		"-in1", implode(" ", $in_for),
		"-in2", implode(" ", $in_rev),
		"-out1 {$fastq_trimmed1}",
		"-out2 {$fastq_trimmed2}",
		"-a1", $sys["adapter1_p5"],
		"-a2", $sys["adapter2_p7"],
		"-qc", $qc_fastq,
		"-threads", bound($threads, 1, 6)
	);
	$parser->exec(get_path("ngs-bits")."SeqPurge", implode(" ", $seqpurge_params), true);
}
else
{
	$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 ".implode(" ", $in_for)." -out $qc_fastq", true);

	$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");

	$skewer_params = array(
		"-x", $sys["adapter1_p5"],
		"-y", $sys["adapter2_p7"],
		"-m any",
		"--threads", $threads,
		"--stdout",
		"--end-quality 15",
		"-"
	);

	$pipline = array();
	$pipeline[] = array("zcat", implode(" ", $in_for));
	$pipeline[] = array(get_path("skewer"), implode(" ", $skewer_params));
	$pipeline[] = array("gzip", "-1 > {$fastq_trimmed1}");
	$parser->execPipeline($pipeline, "skewer");
}

//run bismark
$bismark_tmp = $parser->tempFolder();
$tmp_out_folder = "{$out_folder}/bismark";
if (!is_dir($tmp_out_folder))
{
	mkdir($tmp_out_folder);
}

//bismark spawns two bowtie processes, use approx. half of specified threads
$threads_parg = max(1, intdiv($threads, 2));

//bismark alignment
$genome = get_path("data_folder")."/genomes/bismark/".$sys['build'];
$bismark_params = [
	"--samtools_path", dirname(get_path("samtools")),
	"--bowtie2 --path_to_bowtie ".dirname(get_path("bowtie2")),
	"--rg_tag --rg_id SAMPLE --rg_sample $name",
	"--basename $name",
	"-N 1",
	"-o $tmp_out_folder",
	"--temp_dir $bismark_tmp",
	"--unmapped",
	"--ambiguous"
];
if ($threads_parg > 1)
{
	$bismark_params[] = "-p $threads_parg";
}
$bismark_params[] = $genome;
if ($paired)
{
	$bismark_params[] = "-1 $fastq_trimmed1 -2 $fastq_trimmed2";
}
else
{
	$bismark_params[] = "$fastq_trimmed1";
}

$parser->exec(get_path("bismark"), implode(" ", $bismark_params), true);

//bismark output files
$bismark_suffix = $paired ? "_pe" : "";
$bismark_bam = $tmp_out_folder."/".$name.$bismark_suffix.".bam";

//bismark deduplication
$args_dedup = [
	"--output_dir", $tmp_out_folder,
	"--bam",
	"--samtools_path", dirname(get_path("samtools")),
];
if ($paired)
{
	$args_dedup[] = "--paired";
}
else
{
	$args_dedup[] = "--single";
}
$args_dedup[] = $bismark_bam;
$parser->exec(dirname(get_path("bismark"))."/deduplicate_bismark", implode(" ", $args_dedup), true);
$bismark_bam_dedup = $tmp_out_folder."/".$name.$bismark_suffix.".deduplicated.bam";

//bismark methylation extractor
$args_extractor = [
	"--samtools_path", dirname(get_path("samtools")),
	"--gzip",
	"--output", $tmp_out_folder,
	"--parallel", $threads,
	"--bedGraph"
];
if ($paired)
{
	$args_extractor[] = "--paired-end --no_overlap";
}
else
{
	$args_extractor[] = "--single-end";
}
$args_extractor[] = $bismark_bam;
$parser->exec(dirname(get_path("bismark"))."/bismark_methylation_extractor", implode(" ", $args_extractor), true);

//bam2nuc
//TODO disabled due to long runtime
//$parser->exec(dirname(get_path("bismark"))."/bam2nuc", implode(" ", [
//	"--samtools_path", dirname(get_path("samtools")),
//	"--dir", $tmp_out_folder,
//	"--genome_folder", $genome,
//	$bismark_bam
//]), true);

//bismark report
$pe = $paired ? "PE" : "SE";
$parser->exec(dirname(get_path("bismark"))."/bismark2report", "--dir $tmp_out_folder --alignment_report {$tmp_out_folder}/{$name}_{$pe}_report.txt", true);
$parser->copyFile($tmp_out_folder."/".$name."_{$pe}_report.html", "{$out_folder}/{$name}_bismark_report.html");

//TODO copy bedgraph and coverage files
//$parser->moveFile($tmp_out_folder."/".$name."_pe.bismark.cov.gz", $out_folder."/".$name."_CpG.tsv.gz");
//$parser->moveFile($tmp_out_folder."/".$name."_pe.bedGraph.gz", $out_folder."/".$name."_CpG.bedGraph.gz");

//TODO handle unmapped reads for accurate mapping statistics?
//$tmp_out_folder."/".$name."_unmapped_reads_1.fq.gz"
//$tmp_out_folder."/".$name."_unmapped_reads_2.fq.gz"

//mark duplicates, sort and index BAM
$out = "{$out_folder}/{$name}.bam";
$pipeline = [];
$pipeline[] = [get_path("samtools")." view", "-h $bismark_bam"];
$pipeline[] = [get_path("samblaster"), ""];
$tmp_for_sorting = $parser->tempFile();
$pipeline[] = [get_path("samtools")." sort", "-T {$tmp_for_sorting} -@ {$threads} -m 1G -o $out"];
$parser->execPipeline($pipeline, "sort and mark duplicates");
$parser->indexBam($out, $threads);

//run mapping QC
$qc_map = "{$out_folder}/{$name}_stats_map.qcML";
$params = [
	"-in ".$out,
	"-out ".$qc_map,
	"-ref ".genome_fasta($sys['build'])
];
if ($sys['target_file'] === "" || $sys['type'] === "WGS")
{
	$params[] = "-wgs";
}
else
{
	$params[] = "-roi ".$sys['target_file'];
}
if(!in_array($sys['build'], [ "hg19", "hg38", "GRCh37", "GRCh38" ]))
{
	$params[] = "-no_cont";
}
$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $params), true);
$parser->execTool("NGS/db_import_qc.php", "-id $name -files $qc_fastq $qc_map -force");
?>