<?php

/**
	@page mapping_subjunc
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("mapping_subjunc", "Alignment of RNA-seq FASTQ files to a reference genome using subjunc.");
//mandatory arguments
$parser->addInfile("in1", "Input forward reads in FASTQ(.GZ) format.", false);
$parser->addOutfile("out", "Output BAM file name.", false);

//optional arguments
$parser->addInfile("in2",  "Input reverse reads in FASTQ(.GZ) format for paired-end alignment.", true);
$parser->addString("gtf",  "Annotation", true, get_path("data_folder")."dbs/gene_annotations/GRCh37.gtf");
$parser->addString("genome", "Genome index", true, get_path("data_folder")."genomes/subread/GRCh37");

$parser->addInt("threads", "Number of parallel threads", true, 8);

extract($parser->parse($argv));

//path to output BAM, and also prefix for additional files
$subjunc_out = $parser->tempFolder()."/aligned";

//build command
$arguments = array(
	"-i $genome",
	"-r $in1",
	"-a $gtf",
	"-T $threads",
	"--allJunctions",
	"-o $subjunc_out"
	);
if(isset($in2))
{
  $arguments[] = "-R $in2";
}

//add a read group line to the final BAM file
$basename = basename($out, ".bam");
$arguments[] = "--rg-id {$basename}";

list($gs, $ps) = explode("_", $basename."_99");
$group_props = array("SM:{$gs}",
	"LB:{$gs}_{$ps}",
	"CN:medical_genetics_tuebingen",
	"DT:".date("c"),
	"PL:ILLUMINA");

$arguments[] = "--rg ".implode(" --rg ", $group_props);

//execute mapping with subjunc
$parser->exec(get_path("subread"), implode(" ", $arguments), true);

//pipeline to mark duplicates and sort
$pipeline = array();

//duplicate flagging with samblaster
$pipeline[] = array(get_path("samtools"), "view -h {$subjunc_out}");
$pipeline[] = array(get_path("samblaster"), "");

//convert SAM to BAM (uncompressed) with samtools
$pipeline[] = array(get_path("samtools"), "view -hu -");

//sort BAM by coordinates
$tmp_for_sorting = $parser->tempFile();
$pipeline[] = array(get_path("samtools"), "sort -T $tmp_for_sorting -m 1G -@ ".min($threads, 4)." -o $out -", true);

//execute (samtools BAM to SAM -> samblaster -> samtools SAM to BAM -> samtools sort)
$parser->execPipeline($pipeline, "mapping");

//keep downstream analysis files
$prefix = substr($out, 0, -4);
$parser->moveFile("{$subjunc_out}.junction.bed", "{$prefix}_junction.bed");
$parser->moveFile("{$subjunc_out}.breakpoints.txt", "{$prefix}_breakpoints.txt");
$parser->moveFile("{$subjunc_out}.indel", "{$prefix}_indel.vcf");

//create BAM index file
$parser->indexBam($out, $threads);

?>