<?php
/** 
	@page mvh_run_qc
*/

require_once("mvh_functions.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("mvh_run_qc", "Run the MVH QC pipeline for the given samples.");
$parser->addInfile("samples", "File with one sample per line.", false);
$parser->addString("out", "TSV file with QC values for the samples.", false);
$parser->addInt("threads", "The maximum number of threads used.", true, 10);
extract($parser->parse($argv));

//init
$db = DB::getInstance("NGSD");
$mvh_folder = get_path("mvh_folder");

//output file does not exist => init with header
if (!is_file($out))
{
	$out_headers = ["#sample", "is_tumor", "sys_name", "sys_type", "GRZ_base_quality_threshold", "percent_bases_above_quality_threshold", "GRZ_depth_cutoff", "sample_mean_depth", "GRZ_min_coverage_cutoff", "percent_roi_above_cutoff"];
	file_put_contents($out, implode("\t", $out_headers)."\n");
}


foreach(file($samples) as $line)
{
	$time_start = microtime(true);
	if (trim($line) == "" || $line[0] == "#") continue;
	
	$ps_name = trim($line); 
	
	//determine processing system
	$system = "";
	$sys = load_system($system, $ps_name);
	$is_wes = $sys['type']=="WES";
	$is_wgs = $sys['type']=="WGS";
	$is_lrgs = $sys['type']=="lrGS";
	$build = $sys['build'];
	$info = get_processed_sample_info($db, $ps_name, false);
	$sys_name = $info['sys_name_short'];

	print "{$ps_name} (".$sys['type']."):\n";
	
	//determine ROI if not WGS/lrGS
	$roi = "";
	if ($sys_name=="twistCustomExomeV2" || $sys_name=="twistCustomExomeV2Covaris") $roi = "{$mvh_folder}/rois/twist_exome_core_plus_refseq.bed";
	if ($roi=="" && !$is_wgs && !$is_lrgs) trigger_error("Could not determine target region for sample '{$ps_name}' with processing system '{$sys_name}'!", E_USER_ERROR);
		
	//generate FASTQs
	$tmp_folder_base = $parser->tempFolder("mvh_run_qc_".$ps_name."_");
	print "  generating FASTQ files in folder {$tmp_folder_base} ...\n";
	$bam = $info['ps_bam'];
	$fq1 = $tmp_folder_base."/{$ps_name}_R1.fastq.gz";
	$fq2 = $tmp_folder_base."/{$ps_name}_R2.fastq.gz";
	$add_args = "";
	if (!$is_lrgs) $add_args = "-out2 {$fq2}";
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$bam} -out1 {$fq1} {$add_args}", [$bam], [$tmp_folder_base]);
	
	
	
	//run QC pipeline
	$is_tumor = $info['is_tumor'];
	$is_somatic = $is_tumor || !$is_wgs;
	$qc_folder = $tmp_folder_base."/qc/";
	exec2("mkdir -p {$qc_folder}");
	$grz_qc = run_qc_pipeline($ps_name, $bam, $fq1, $fq2, $roi, $is_tumor, $parser, $info['sys_type'], $is_somatic, $qc_folder, $is_lrgs, "DUMMY_ID", $threads);
	$minQual  = $grz_qc["qualityThreshold"];
	$percQual = number_format($grz_qc["percentBasesAboveQualityThreshold"], 2);
	$minDepth = $grz_qc["meanDepthOfCoverageRequired"];
	$meanDepth = (float)($grz_qc["meanDepthOfCoverage"]);
	$minCov = $grz_qc["minCoverage"];
	$regionsAboveMin = number_format(100.0*$grz_qc["targetedRegionsAboveMinCoverage"], 2);
	
	//output to file
	$out_line = implode("\t", [$ps_name, $is_tumor ? "yes" : "no", $sys_name, $info['sys_type'], $minQual, $percQual, $minDepth, $meanDepth, $minCov, $regionsAboveMin])."\n";
	file_put_contents($out, $out_line, FILE_APPEND);
	
	//output to command line
	print "  overall runtime: ".time_readable(microtime(true)-$time_start)."\n";;
	print "  base quality cutoff: {$minQual} - percentage bases above cutoff: {$percQual}\n";
	print "  min depth: {$minDepth} - mean depth: {$meanDepth} - mean depth plus 5%: ".number_format($meanDepth*1.05, 2)."\n";
	print "  min coverage: {$minCov} - percentage target regions above min coverage: {$regionsAboveMin}\n";
	
	//remove temp files when finished
	exec2("rm -rf {$tmp_folder_base}");
}

?>