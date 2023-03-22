<?php 
/** 
	@page vc_cnvkit
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_cnvkit", "Copy-number calling with CNVkit.");

$parser->addInfile("bam",  "Input BAM file.", true);
$parser->addString("output_folder", "Data folder to store the analysis.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
//bam (references)

//optional
$parser->addInfile("target_region",  "Enrichment targets BED file.", true);
$parser->addInfile("exclude_region",  "BED file containing unmappable/variable/poorly sequenced regions which should be excluded.", true);
$parser->addInt("binsize", "", true, 267);
$parser->addInt("min_binsize", "", true, 267);
$parser->addInt("min_gapsize", "", true, -1);
//center/method/ploidy for cnvkit call?
//params for plot
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//init
if (!is_dir($output_folder))
{
	trigger_error("Output folder '" + $output_folder + "' doesn't exist!");
}
$basename = $output_folder."/".$name;

$cnvkit = dirname(get_path("python3"))."/cnvkit.py";
$genome = genome_fasta($build);

//create targets
$targets = $basename."_cnvs_targets.bed";
$parser->exec($cnvkit, "target {$target_region} --avg-size {$binsize} --split --output {$targets}");

//access
$access_bed_file = $basename."_cnvs_access.bed";
$parser->exec($cnvkit, "access {$genome} --exclude {$exclude_region} --min-gap-size {$min_gapsize} --output {$access_bed_file}");

//create anti-targets
$anti_targets = $basename."_cnvs_anti_targets.bed";
$parser->exec($cnvkit, "antitarget {$target_region} --access {$access_bed_file} --avg-size {$binsize} --min-size {$min_binsize} --output {$anti_targets}");

//calculate target coverage
$target_coverage = $basename."_cnvs_target_coverage.cnn";
$parser->exec($cnvkit, "coverage {$bam} {$targets} --processes {$threads} --output {$target_coverage}");

//calculate anti-target coverage
$anti_target_coverage = $basename."_cnvs_anti_target_coverage.cnn";
$parser->exec($cnvkit, "coverage {$bam} {$anti_targets} --processes {$threads} --output {$anti_target_coverage}");

//create reference
//TODO: collect reference samples
$cohort_bams = array();
$copy_number_reference = $basename."_cnvs_reference.cnn";
$parser->exec($cnvkit, "reference --output {$copy_number_reference}"); //TODO

//fix coverage
$copy_number_ratio_table = $basename."_cnvs_ratio_table.cnr";
$parser->exec($cnvkit, "fix {$target_coverage} {$anti_target_coverage} {$copy_number_reference} --processes {$threads} --output {$copy_number_ratio_table}");

//calculate discrete copy-numbers
$copy_number_raw_file = $basename."_cnvs_raw.cns";
$method = ""; //TODO: determine segmentation method
$parser->exec($cnvkit, "segment {$copy_number_ratio_table} --method {$method} --processes {$threads} --output {$copy_number_raw_file}");

$copy_number_call_file = $basename."_raw_cnvs_cnvkit.cns";
$parser->exec($cnvkit, "call {$copy_number_raw_file} --center {$center} --method {$method} --ploidy {$ploidy} --processes {$threads} --output {$copy_number_call_file}");

//create plots
$copy_number_plot = $basename."_cnvs_plot.pdf";
$region = ""; //TODO: implement
$color = "red";
$y_max = 10;
$y_min = 0;
$parser->exec($cnvkit, "scatter {$copy_number_raw_file} --segment {$copy_number_raw_file} --y-max {$y_max} --y-min {$y_min} --output {$copy_number_plot}");

//determine gender
$gender_file = $basename."_cnv_gender.txt";
$parser->exec($cnvkit, "sex {$copy_number_reference} {$copy_number_ratio_table} {} --output {$gender_file}");



//create pdf





?>
