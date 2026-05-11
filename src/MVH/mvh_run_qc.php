<?php
/** 
	@page mvh_run_qc
*/

require_once("mvh_functions.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("mvh_run_qc", "Run the QC pipeline for the given samples.");
$parser->addInfile("samples", "File with one sample per line.", false);
$parser->addString("out", "TSV file with QC values for the samples.", false);
extract($parser->parse($argv));

$out_headers = ["#sample", "is_tumor", "sys_name", "sys_type", "GRZ_quality_threshold", "GRZ_percent_above_qality_threshold", "GRZ_mean_depth", "GRZ_mean_depth_plus_5percent", "GRZ_min_coverage_regions", "GRZ_percent_regions_above_region_threshold"];

if (! is_file($out)) file_put_contents($out, implode("\t", $out_headers)."\n");


foreach(file($samples) as $line)
{
	$time_start = microtime(true);
	if (trim($line) == "" || $line[0] == "#") continue;
	
	$ps_name = trim($line); 
	//determine processing system
	$sys = load_system($system, $ps_name);
	$is_wes = $sys['type']=="WES";
	$is_wgs = $sys['type']=="WGS";
	$roi = trim($sys['target_file']);
	$build = $sys['build'];
	
	$db = DB::getInstance("NGSD");
	$info = get_processed_sample_info($db, $ps_name, false);
	$sys_name = $info['sys_name_short'];
	$sys_type = $info['sys_type'];
	
	$seq_mode = $sys_type;
	$is_lrgs = $sys_type=="lrGS";
	
	
	if ($is_lrgs) trigger_error("Long read samples currently not supported (implementation ofthe BAM correction necessary?): $ps_name", E_USER_ERROR);
	
	$is_tumor = $info['is_tumor'];
	$is_somatic = $is_tumor;

	$tmp_folder_base = $parser->tempFolder("mvh_run_qc_".$ps_name."_");
	$qc_folder = $tmp_folder_base."/qc/";
	exec2("mkdir -p {$qc_folder}");
	
	
	$bam = $info['ps_bam'];
	$fq1 = $tmp_folder_base."/{$ps_name}_R1.fastq.gz";
	$fq2 = $tmp_folder_base."/{$ps_name}_R2.fastq.gz";
	
	print "  generating FASTQ files for sample {$ps_name}...\n";
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$bam} -out1 {$fq1} -out2 {$fq2}", [$bam], ["{$tmp_folder_base}/"]);
	
	$patient_id = "DUMMY_ID";
	
	$grz_qc = run_qc_pipeline($ps_name, $bam, $fq1, $fq2, $roi, $is_tumor);
	print "QC: $ps_name - time taken: ".time_readable(microtime(true)-$time_start)."\n";;
	print_grz_qc($grz_qc);
	
	$minQual  = $grz_qc["qualityThreshold"];
	$percQual = $grz_qc["percentBasesAboveQualityThreshold"];
	$meanDepth = (float)($grz_qc["meanDepthOfCoverage"]);
	$minCov = $grz_qc["minCoverage"];
	$regionsAboveMin = number_format($grz_qc["targetedRegionsAboveMinCoverage"],2);
	
	$out_line = implode("\t", [$ps_name, $is_tumor ? "yes" : "no", $sys_name, $sys_type, $minQual, $percQual, $meanDepth, number_format($meanDepth*1.05,2), $minCov, $regionsAboveMin])."\n";
	
	file_put_contents($out, $out_line, FILE_APPEND);
	
	//remove temp files when finished
	exec2("rm -rf {$tmp_folder_base}");
	
}

?>