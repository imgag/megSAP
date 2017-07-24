<?php

/*
	@page somatic_pair_rna
 
	@todo update old fastq handling
 */

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";
require_once($basedir."Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("somatic_rna", "Analysis tumor normal RNA samples.");
$parser->addString("p_folder","Folder containing sample subfolders with fastqs (Sample_GSXYZ).",false);
$parser->addString("t_id", "Tumor sample processing-ID (e.g. GSxyz_01). There should be a folder 'Sample_tsid' that contains fastq-files marked with the id within the p_folder.", false);
$parser->addString("o_folder", "Output folder.", false);
//optional
$parser->addString("n_id", "Normal sample processing-ID (e.g. GSxyz_01). There should be a folder 'Sample_nsid' that contains fastq-files marked with the id within the p_folder.", true, "na");
$parser->addString("t_folder", "Folder where original tumor fastq files can be found.", true, "na");
$parser->addString("n_folder", "Folder where original normal fastq files can be found.", true, "na");
$parser->addInfile("t_sys",  "Tumor processing system INI file (determined from 't_id' by default).", true);
$parser->addInfile("n_sys",  "Reference processing system INI file (determined from 'n_id' by default).", true);
$steps_all = array("ma","fu","vc","an","db");
$parser->addString("steps", "Comma-separated list of processing steps to perform. (".implode(",",$steps_all).")", true, implode(",", array_slice($steps_all,1)));
$parser->addFlag("nsc", "No sample correlation check.");
extract($parser->parse($argv));

$parser->log("Pipeline revision: ".repository_revision(true));

// (0) preparations
// (0a) determine steps to perform
// (0b) check if tumor-normal pair or not
$tumor_only = false;
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!",E_USER_ERROR);
}
if($n_id=="na")
{
	$tumor_only = true;
	$available_steps = array("fastq","ma","db");
	$steps = array_intersect($available_steps,$steps);
}
//$steps_perf = array_slice($steps_all, $start_index, $end_index-$start_index+1);
// (0c) check in_folder
$t_folder = $p_folder."/Sample_".$t_id."/";
$n_folder = $p_folder."/Sample_".$n_id."/";
// (00) check out put folder
$o_folder = $o_folder."/";
// (0e) get tum and ref systems
$o_folder = $o_folder."/";
$o_folder_tum = $t_folder;
$o_folder_ref = $n_folder;
if (!file_exists($o_folder))	trigger_error("Output-folder '$o_folder' does not exist.", E_USER_ERROR);
if (!file_exists($o_folder_tum))	mkdir($o_folder_tum);
if (!$tumor_only && !file_exists($o_folder_ref))	mkdir($o_folder_ref);
if (!file_exists($t_folder))	trigger_error("Tumor-folder '$t_folder' does not exist.", E_USER_ERROR);
if (!$tumor_only && !file_exists($n_folder))	trigger_error("Reference-folder '$n_folder' does not exist.", E_USER_ERROR);

$t_sys_ini = load_system($t_sys, $t_id);
// copy tumor ini file to folder that is accessible from all servers
$tmp = $t_sys;
$t_sys = $t_folder."system.ini";
if(!copy($tmp,$t_sys))	trigger_error("Could not copy ini file from $tmp to $t_sys.",E_USER_ERROR);
if(!$tumor_only)
{
	$n_sys_ini = load_system($n_sys, $n_id);
	// copy normal ini file to folder that is accessible from all servers
	$tmp = $n_sys;
	$n_sys = $n_folder."system.ini";
	if(!copy($tmp,$n_sys))	trigger_error("Could not copy ini file from $tmp to $n_sys.",E_USER_ERROR);
}
if(!$tumor_only && $t_sys_ini['name_short'] != $n_sys_ini['name_short']) trigger_error ("System tumor '".$t_sys_ini['name_short']."' and system reference '".$n_sys_ini['name_short']."' are different!", E_USER_WARNING);
if(!$tumor_only && $t_sys_ini['build'] != $n_sys_ini['build']) trigger_error ("System tumor '".$t_sys_ini['build']."' and system reference '".$n_sys_ini['build']."' do have different builds!", E_USER_ERROR);

// (0) set RNA fastq files
$t_forward = $o_folder_tum.$t_id."_*_R1_001.fastq.gz";
$t_reverse = $o_folder_tum.$t_id."_*_R2_001.fastq.gz";
if(is_valid_processingid($t_id))
{
	list($s_id,$p_id) = explode("_",$t_id);
	$t_forward = $o_folder_tum.$s_id."_*_R1_001.fastq.gz";
	$t_reverse = $o_folder_tum.$s_id."_*_R2_001.fastq.gz";
}
$n_forward = $o_folder_ref.$n_id."_*_R1_001.fastq.gz";
$n_reverse = $o_folder_ref.$n_id."_*_R2_001.fastq.gz";
if(is_valid_processingid($n_id))
{
	list($s_id,$p_id) = explode("_",$n_id);
	$n_forward = $o_folder_ref.$s_id."_*_R1_001.fastq.gz";
	$n_reverse = $o_folder_ref.$s_id."_*_R2_001.fastq.gz";
}

// (1) map reference and tumor sample
$tum_bam = $o_folder_tum.$t_id.".bam";
$ref_bam = $o_folder_ref.$n_id.".bam";
$tum_counts = $o_folder_tum.$t_id."_counts_fpkm.tsv";
$ref_counts = $o_folder_ref.$n_id."_counts_fpkm.tsv";
if(in_array("ma", $steps))
{	
	//map tumor and normal in high-mem-queue
	$args = "-steps ma,rc,fu,an";
	$working_directory = realpath($p_folder);
	$commands = array("php ".$basedir."Pipelines/analyze_rna.php -in_for $t_forward -in_rev $t_reverse -system $t_sys -out_folder ".$o_folder_tum." -out_name $t_id $args --log ".$o_folder_tum."analyze_".date('YmdHis',mktime()).".log");
	if(!$tumor_only)	$commands[] = "php ".$basedir."Pipelines/analyze_rna.php -in_for $n_forward -in_rev $n_reverse -system $n_sys -out_folder ".$o_folder_ref." -out_name $n_id $args --log ".$o_folder_ref."analyze_".date('YmdHis',mktime()).".log";
	$parser->jobsSubmit($commands, $working_directory, get_path("queues_high_mem"), true);
}

// only for tumor-normal pairs
if(!$tumor_only)
{
	// calculate counts tumor, normal and somatic fold change
	$som_counts = $o_folder.$t_id."-".$n_id."_counts_fpkm.tsv";
	if(!$tumor_only)	$parser->execTool("NGS/rc_compare.php", "-in1 $tum_counts -in2 $ref_counts -out $som_counts");

	// (2) check that samples are related
	if(!$nsc)
	{
		$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in1 $tum_bam -in2 $ref_bam -bam -max_snps 4000", true);
		$parts = explode(":", $output[0][1]);
		if ($parts[1]<0.8)
		{
			trigger_error("The genotype correlation of samples $tum_bam and $ref_bam is ".$parts[1]."; it should be above 0.8!", E_USER_ERROR);
		}
	}

	$tum_fu = $o_folder_tum.$t_id."_var_fusions.tsv";
	$nor_fu = $o_folder_ref.$n_id."_var_fusions.tsv";
	$som_fu = $o_folder.$t_id."-".$n_id."_var_fusions.tsv";
	if(in_array("fu", $steps))
	{	
		$fusions1 = Matrix::fromTSV($tum_fu);	//tumor
		$idx_tleft = $fusions1->getColumnIndex("LeftBreakpoint");
		$idx_tright = $fusions1->getColumnIndex("RightBreakpoint");
		$fusions2 = Matrix::fromTSV($nor_fu);	//normal
		$idx_nleft = $fusions2->getColumnIndex("LeftBreakpoint");
		$idx_nright = $fusions2->getColumnIndex("RightBreakpoint");
		
		$fusions_somatic = new Matrix();
		$fusions_somatic->setHeaders($fusions1->getHeaders());
		for($i=0;$i<$fusions1->rows();++$i)
		{
			$somatic = true;

			$r_tum = $fusions1->getRow($i);
			
			for($j=0;$j<$fusions2->rows();++$j)
			{
				$r_nor = $fusions2->getRow($j);
				
				if($r_tum[$idx_tleft]==$r_nor[$idx_nleft] && $r_tum[$idx_tright]==$r_nor[$idx_nright])	$somatic = false;
			}
			
			if($somatic)	$fusions_somatic->addRow($r_tum);
		}
		$fusions_somatic->toTSV($som_fu);
	}

	// (3) run strelka
	$som_v = $o_folder.$t_id."-".$n_id."_var.vcf.gz";
	$som_vcf = $o_folder.$t_id."-".$n_id."_var_annotated.vcf.gz";
	if(in_array("vc", $steps))
	{
		// (3a) variant calling
		$args = "";
		if (!$t_sys_ini['shotgun']) $args .= "-amplicon ";
		$parser->execTool("NGS/vc_strelka.php", "-t_bam $tum_bam -n_bam $ref_bam -out $som_v $args");
	}


	// (4) annotation
	$som_gsvar = $o_folder.$t_id."-".$n_id.".GSvar";
	if(in_array("an", $steps))
	{
		// (4a) annotate vcf
		$tmp1 = $t_id;
		$tmp2 = $n_id;
		if(is_valid_processingid($t_id))	list($tmp1,) = explode("_", $t_id);	
		if(is_valid_processingid($n_id))	list($tmp2,) = explode("_", $n_id);	
		$parser->execTool("Pipelines/annotate.php", "-out_name ".basename($som_gsvar, ".GSvar")." -out_folder ".dirname($som_gsvar)." -system ".$t_sys." -vcf $som_v -t_col $tmp1 -n_col $tmp2");
		
		// (4b)	convert to GSvar
		$extra = "-t_col $t_id ";
		if($n_id!="na")	 $extra .= "-n_col $n_id";
		$parser->execTool("NGS/vcf2gsvar_somatic.php", "-in $som_vcf -out $som_gsvar $extra");
		$parser->execTool("NGS/an_dbNFSPgene.php", "-in $som_gsvar -out $som_gsvar");
		
		// (4c) Annotate somatic NGSD-data
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode somatic", true);
		// (4d) annotate frequencies and depths
		//global settings for annotation
		$vaf_options = " -depth";
		// re-annotate depth and variant frequency tumor (to make sure we used the same version of VariantAnnotateFrequency both for tumor and normal tissue)
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $som_gsvar -bam $tum_bam -out $som_gsvar -name dna_tum $vaf_options", true);
		// re-annotate depth and variant frequency reference
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $som_gsvar -bam $ref_bam -out $som_gsvar -name dna_ref $vaf_options", true);
		// filter somatic variant list for somatic variants
		//TODO $parser->execTool("NGS/filter_tsv.php", "-in $som_gsvar -out $som_gsvar -type somatic -roi ".$n_sys_ini['target_region']);
	}
}

//	(5)	import statistics to NGSD
if (in_array("db", $steps))
{
	//	(4a)	tumor sample
	if(is_valid_processingid($t_id))
	{
		//$parser->execTool("NGS/db_check_gender.php", "-in $t_bam -pid $t_id");

		//import QC data tumor
		$log_db  = $t_folder."/".$t_id."_log4_db.log";
		$qc_fastq  = $t_folder."/".$t_id."_stats_fastq.qcML";
		$qc_map  = $t_folder."/".$t_id."_stats_map.qcML";	//MappingQC does currently not support RNA
		$parser->execTool("NGS/db_import_qc.php", "-id $t_id -files $qc_fastq -force -min_depth 0 --log $log_db");

		//update last_analysis date
		updateLastAnalysisDate($t_id, $tum_bam);
		if(!isTumor($t_id))	trigger_error("Tumor $t_id is not flagged as tumor in NGSD",E_USER_WARNING);
	}
	else	trigger_error("No DB import since no valid processing ID (T:".$t_id.")",E_USER_WARNING);

	//	(4b)	normal sample
	if(!$tumor_only)
	{
		if(is_valid_processingid($n_id))
		{
			//$parser->execTool("NGS/db_check_gender.php", "-in $n_bam -pid $n_id");

			//import QC data normal
			$log_db  = $n_folder."/".$n_id."_log4_db.log";
			$qc_fastq  = $n_folder."/".$n_id."_stats_fastq.qcML";
			$qc_map  = $n_folder."/".$n_id."_stats_map.qcML";	//MappingQC does currently not support RNA
			$parser->execTool("NGS/db_import_qc.php","-id $n_id -files $qc_fastq -force -min_depth 0 --log $log_db");

			//update last_analysis date
			updateLastAnalysisDate($n_id, $ref_bam);
			
			// update normal entry for tumor
			if(updateNormalSample($t_id, $n_id))	trigger_error("Updated normal sample ($n_id) for tumor ($t_id) within NGSD.",E_USER_NOTICE);

			if(isTumor($n_id))	trigger_error("Normal $n_id is flagged as tumor in NGSD",E_USER_WARNING);
		}
		else	trigger_error("No DB import since no valid processing ID (N:".$n_id.")",E_USER_WARNING);
	}
}
?>