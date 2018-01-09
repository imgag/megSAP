<?php

/**
	@page somatic_dna
	@todo recalculate QC-values at the end of the pipeline
	@todo consider replacing analyze.php by own combinations of tools
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_dna", "Differential analysis of sequenced tumor and normal DNA samples.");
$parser->addString("p_folder", "Folder that contains Sample folders.", false);
$parser->addString("t_id", "Tumor sample processing-ID (e.g. GSxyz_01).", false);
$parser->addString("n_id", "Normal sample processing-ID (e.g. GSxyz_01). To process a tumor samples solely use 'na'.", false);
$parser->addString("o_folder", "Output folder.", false);
$steps_all = array("ma", "vc", "an", "ci", "db");
$parser->addString("steps", "Comma-separated list of processing steps to perform. Available are: ".implode(",", $steps_all), true, "ma,vc,an,db");
$parser->addString("filter_set","Filter set to use. Only if annotation step is selected. Multiple filters can be comma separated.",true,"synonymous,not-coding-splicing");
// optional
$parser->addFloat("min_af", "Allele frequency detection limit.", true, 0.05);
$parser->addInfile("t_sys",  "Tumor processing system INI file (determined from the NGSD using 't_id' by default).", true);
$parser->addInfile("n_sys",  "Reference processing system INI file (determined from the NGSD using 'n_id' by default).", true);
$parser->addFlag("abra", "Use Abra for combined indel realignment.");
$parser->addFlag("smt", "Skip mapping tumor (only in pair mode).");
$parser->addFlag("smn", "Skip mapping normal (only in pair mode).");
$parser->addFlag("nsc", "Skip sample correlation check (human only, only in pair mode).");
$parser->addFlag("keep_all_variants_strelka", "Reduce number of variants to PASS by strelka.");
$parser->addFlag("reduce_variants_filter", "Reduce number of variants to PASS by own filter tool.");
$parser->addFlag("freebayes", "Use freebayes for variant calling (default: strelka).");
$parser->addFloat("contamination", "Indicates fraction of tumor cells in normal sample.", true, 0);
$parser->addFlag("no_softclip", "Skip soft-clipping of overlapping reads. NOTE: This may increase rate of false positive variants.", true);
$parser->addInfile("promoter","Bed file containing promoter regions. Will be used for filter column of vcf.",true);
$parser->addEnum("clip", "Soft-clip overlapping read pairs.", true, array("sc","mfb","mfm","mfr"),"sc");
$parser->addFlag("strelka1","Use strelka1 for variant calling.",true);
$parser->addFlag("add_vc_folder","Add folder containing variant calling results from variant caller.",true);
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));

// (0) preparations
// check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

// tumor only or tumor normal
$single_sample = ($n_id === "na");
if($single_sample)	trigger_error("Single sample mode.",E_USER_NOTICE);
else	trigger_error("Paired sample mode.",E_USER_NOTICE);

// check and generate output folders
$o_folder = $o_folder."/";
$t_folder = $p_folder."/Sample_".$t_id."/";
if (!file_exists($t_folder))	trigger_error("Tumor-folder '$t_folder' does not exist.", E_USER_ERROR);
$t_sys_ini = load_system($t_sys, $t_id);
$n_folder = $p_folder."/Sample_".$n_id."/";
if (!$single_sample && !file_exists($n_folder))	trigger_error("Reference-folder '$n_folder' does not exist.", E_USER_ERROR);
// get ref system and do some basic checks
if(!$single_sample)
{
	$n_sys_ini = load_system(($n_sys), $n_id);
	if(empty($t_sys_ini['target_file']) || empty($n_sys_ini['target_file'])) trigger_error ("System tumor or system normal does not contain a target file.", E_USER_WARNING);
	if($t_sys_ini['name_short'] != $n_sys_ini['name_short']) trigger_error ("System tumor '".$t_sys_ini['name_short']."' and system reference '".$n_sys_ini['name_short']."' are different!", E_USER_WARNING);
	if($t_sys_ini['build'] != $n_sys_ini['build']) trigger_error ("System tumor '".$t_sys_ini['build']."' and system reference '".$n_sys_ini['build']."' do have different builds!", E_USER_ERROR);
}

//amplicon mode
$amplicon = false;
if(!$t_sys_ini['shotgun'])
{
	$amplicon = true;
	trigger_error("Amplicon mode.",E_USER_NOTICE);
}

// mapping
$t_bam = $t_folder.$t_id.".bam";
$n_bam = $n_folder.$n_id.".bam";
if (in_array("ma", $steps))
{
	// run analyze.php on tumor and normal
	$analyze_args = [
		"-steps ma",
		"-no_abra" // disable ABRA realignment due to manta compatibility issues
	];

	// tumor sample
	if (!$smt)
	{
		$analyze_args_tum = [
			"-folder", $t_folder,
			"-name", $t_id,
			"-system", $t_sys,
			"--log", "{$t_folder}analyze_" . date('YmdHis',mktime()) . ".log"
		];
		$parser->execTool("Pipelines/analyze.php", implode(" ", array_merge($analyze_args, $analyze_args_tum)));
	}

	// normal sample
	if (!$single_sample && !$smn)
	{
		$analyze_args_nor = [
			"-folder", $n_folder,
			"-name", $n_id,
			"-system", $n_sys,
			"--log", "{$n_folder}analyze_" . date('YmdHis',mktime()) . ".log"
		];
		$parser->execTool("Pipelines/analyze.php", implode(" ", array_merge($analyze_args, $analyze_args_nor)));
	}

	// overlap clipping, is not done by analyze.php to prevent import of QC data calculated on BamClipOverlap result
	if(!$no_softclip)
	{
		$clip_arg_map = [
			"mfb" => "-overlap_mismatch_baseq",
			"mfm" => "-overlap_mismatch_mapq",
			"mfr" => "-overlap_mismatch_remove",
			"sc" => ""
		];
		$clip_arg = $clip_arg_map[$clip];

		if (!$smt)
		{
			$tmp1_t_bam = $parser->tempFile("_tumor.bam");
			$parser->exec(get_path("ngs-bits")."BamClipOverlap", "-in {$t_bam} -out {$tmp1_t_bam} {$clip_arg}", true);
			$parser->sortBam($tmp1_t_bam, $t_bam, $threads);
			$parser->indexBam($t_bam, $threads);
		}

		if(!$single_sample && !$smn)
		{
			$tmp1_n_bam = $parser->tempFile("_normal.bam");
			$parser->exec(get_path("ngs-bits")."BamClipOverlap", "-in {$n_bam} -out {$tmp1_n_bam} {$clip_arg}", true);
			$parser->sortBam($tmp1_n_bam, $n_bam, $threads);
			$parser->indexBam($n_bam, $threads);
		}
	}

	// combined indel realignment with ABRA
	if ($abra && $t_sys_ini['type']!="WGS" && ($single_sample || $n_sys_ini['type']!="WGS"))
	{
		// intersect both target files
		$tmp_targets = $t_sys_ini['target_file'];
		if(!$single_sample)
		{
			$tmp_targets = $parser->tempFile("_intersected.bed");
			$commands_params = [];
			$commands_params[] = [get_path("ngs-bits")."BedIntersect", "-in ".$t_sys_ini['target_file']." -in2 ".$n_sys_ini['target_file']." -mode intersect"];
			$commands_params[] = [get_path("ngs-bits")."BedMerge", "-out $tmp_targets"];
			$parser->execPipeline($commands_params, "intersect target regions", true);
		}

		// indel realignment with ABRA
		$tmp_in = (!$single_sample?$n_bam." ":"")."$t_bam";
		$tmp1_t_bam = $parser->tempFile("_tumor.bam");
		$tmp1_n_bam = $parser->tempFile("_normal.bam");
		$tmp_out = (!$single_sample?$tmp1_n_bam." ":"")."$tmp1_t_bam";
		$command = "php ".repository_basedir()."/src/NGS/indel_realign_abra.php -in $tmp_in -out $tmp_out -roi $tmp_targets -threads $threads -mer 0.02 -mad 5000 2>&1";
		$working_directory = realpath($p_folder);
		$parser->jobsSubmit(array($command), $working_directory, get_path("queues_high_mem"), true);
		
		// copy realigned files to output folder and overwrite previous bam files
		$parser->moveFile($tmp1_t_bam, $t_bam);
		$parser->indexBam($t_bam, $threads);
		if(!$single_sample)
		{
			$parser->moveFile($tmp1_n_bam, $n_bam);
			$parser->indexBam($n_bam, $threads);
		}
	}	
}

// check if samples are related
if(!$nsc && !$single_sample)
{
	$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in $t_bam $n_bam -mode bam -max_snps 4000", true);
	$correlation = explode("\t", $output[0][1])[3];
	if ($correlation<0.8)
	{
		trigger_error("The genotype correlation of samples $t_bam and $n_bam is {$correlation}; should be above 0.8!", E_USER_ERROR);
	}
}

// variant calling
$som_prefix = $t_id;
$som_v = $o_folder.$t_id."_var.vcf.gz";
if(!$single_sample)	$som_v = $o_folder.$t_id."-".$n_id."_var.vcf.gz";
$som_cnv = $o_folder.$t_id."_var_cnvs.tsv";
if(!$single_sample)	$som_cnv = $o_folder.$t_id."-".$n_id."_cnvs.tsv";
$som_gsvar = $o_folder.$som_prefix.".GSvar";
if(!$single_sample)	$som_qci = $o_folder.$t_id."-".$n_id."_var_qci.vcf.gz";
if(!$single_sample)	$som_sv = $o_folder.$t_id."-".$n_id."_var_structural.vcf.gz";
if(!$single_sample)	$som_svt = $o_folder.$t_id."-".$n_id."_var_structural.tsv";
if(!$single_sample)	$som_si = $o_folder.$t_id."-".$n_id."_var_smallIndels.vcf.gz";
$som_bafs = $o_folder.$t_id."-".$n_id."_bafs.seg";
if (in_array("vc", $steps))
{
	// structural variant calling, should be done before variant calling since strelka uses the smallIndel output
	if(!$single_sample)
	{
		if(!starts_with($t_sys_ini['name_manufacturer'],"ThruPlex TagSeq"))
		{
			if($t_sys_ini['shotgun']==1)
			{
				$par = "";
				if($t_sys_ini['type']=="WES")	$par .= "-exome ";
				if($add_vc_folder)	$par .= "-temp ".dirname($som_v)."/variant_calling ";
				$par .= "-smallIndels $som_si ";
				$parser->execTool("NGS/vc_manta.php", "-t_bam $t_bam -bam $n_bam $par -out $som_sv -build ".$t_sys_ini['build']);
				$parser->execTool("Tools/converter_manta2tsv.php", "-in $som_sv -out $som_svt -tumor_id $t_id");
			}
		}
		else trigger_error("Breakpoint detection deactivated for ThruPlex samples.",E_USER_NOTICE);
	}
	else	trigger_error("Breakpoint detection currently only implemented for tumor normal pairs.",E_USER_NOTICE);
	
	// variant calling
	if($freebayes || $single_sample)	// variant calling using freebayes
	{
		// vc_freebayes uses ngs-bits that can not handel multi-sample vcfs
		$tmp1 = $parser->tempFile();
		$tmp_in = $t_bam.($single_sample?"":" $n_bam");
		$parser->execTool("NGS/vc_freebayes.php","-bam $tmp_in -out $tmp1 -build ".$t_sys_ini['build']." -min_af ".$min_af." -target ".$t_sys_ini['target_file']);

		// fix header
		$tmp2 = $parser->tempFile();
		$s = Matrix::fromTSV($tmp1);
		$tmp_headers = $s->getHeaders();
		$tumor_col = NULL;
		$normal_col = NULL;
		for($i=0;$i<count($tmp_headers);++$i)
		{
			if(contains($tmp_headers[$i],$t_id))
			{
				$tumor_col = $i;
				$tmp_headers[$tumor_col] = $t_id;
			}
			if(!$single_sample && contains($tmp_headers[$i],$n_id))
			{
				$normal_col = $i;
				$tmp_headers[$normal_col] = $n_id;
			}
		}
		$s->setHeaders($tmp_headers);
		//	fix genotype information
		//	fields: GT:DP:RO:QR:AO:QA:GL
		//	is often: .
		//	should be: 1/1:3:0:0:3:103:-9.61333,-0.90309,0
		//	using default: ./.:0:.:.:.:0:.,.,.
		//	@TODO simply skip variants with no information in tumor or normal?
		for($i=0;$i<$s->rows();++$i)
		{
			$row = $s->getRow($i);
			if($row[$tumor_col]==".")
			{
				$s->set($i,$tumor_col,"./.:0:.:.:.:0:.,.,.");
			}
			if(!$single_sample && $row[$normal_col]==".")
			{
				$s->set($i,$normal_col,"./.:0:.:.:.:0:.,.,.");
			}
		}
		// add tumor + normal information
		if($single_sample)	$s->addComment("#PEDIGREE=<Tumor=$t_id>");
		else	$s->addComment("#PEDIGREE=<Tumor=$t_id,Normal=$n_id>");	// add pedigree information for tumor and normal
		$s->toTSV($tmp2);
		// zip and index output file
		$parser->exec("bgzip", "-c $tmp2 > $som_v", false);	// no output logging, because Toolbase::extractVersion() does not return
		$parser->exec("tabix", "-f -p vcf $som_v", false);	// no output logging, because Toolbase::extractVersion() does not return

		// filter for variants that are likely to be germline, included in strelka2					
		$extra = array("-type somatic-lq","-keep");
		$parser->execTool("NGS/filter_vcf.php", "-in $som_v -out $som_v -min_af $min_af  ".implode(" ", $extra));
	}
	else
	{
		$args = array();
		$args[] = "-build ".$t_sys_ini['build'];
		$args[] = "-target ".$t_sys_ini['target_file'];
		if ($keep_all_variants_strelka) $args[] = "-k";
		//if ($amplicon) $args[] = "-amplicon";
		if($add_vc_folder)	$args[] = "-temp ".dirname($som_v)."/variant_calling";
		if(is_file($som_si))	$args[] = "-smallIndels $som_si";
		$parser->execTool("NGS/vc_strelka2.php", "-t_bam $t_bam -n_bam $n_bam -out $som_v ".implode(" ", $args));
		//else 	$parser->execTool("NGS/vc_strelka.php", "-t_bam $t_bam -n_bam $n_bam -out $som_v ".implode(" ", $args));	//combined variant calling using strelka		
	}

	// copy number variant calling
	$tmp_folder = $parser->tempFolder();
	$t_cov = $tmp_folder."/".basename($t_bam,".bam").".cov";
	$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $t_bam -in ".$t_sys_ini['target_file']." -out $t_cov",true);
	if(!$single_sample)	$n_cov = $tmp_folder."/{$n_id}.cov";
	if(!$single_sample)	$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $n_bam -in ".$n_sys_ini['target_file']." -out $n_cov", true);
	//copy coverage file to reference folder (has to be done before CnvHunter call to avoid analyzing the same sample twice)
	if (!$single_sample && is_valid_ref_sample_for_cnv_analysis($n_id))
	{
		//create reference folder if it does not exist
		$ref_folder = get_path("data_folder")."/coverage/".$n_sys_ini['name_short']."/";
		if (!is_dir($ref_folder)) mkdir($ref_folder);			
		//copy file
		$ref_file = $ref_folder.$n_id.".cov";
		$parser->copyFile($n_cov, $ref_file);
		$n_cov = $ref_file;
	}
	$tmp_in = "-cov $t_cov ".($single_sample?"":"-n_cov $n_cov");
	$n = 20;
	if($t_sys_ini['type']=="WES")	$n = 30;
	$parser->execTool("NGS/vc_cnvhunter.php", "$tmp_in -out $som_cnv -system $t_sys -min_corr 0 -seg $t_id -n $n");
	
	// add somatic baf file
	$parser->execTool("NGS/mapping_baf.php", "-in $t_bam -out $som_bafs -system $t_sys");
}

// annotation
$som_gsvar = $o_folder.$t_id.".GSvar";
if(!$single_sample)	$som_gsvar = $o_folder.$t_id."-".$n_id.".GSvar";
$som_unfi = $o_folder.$t_id."_var_annotated_unfiltered.vcf.gz";
if(!$single_sample)	$som_unfi = $o_folder.$t_id."-".$n_id."_var_annotated_unfiltered.vcf.gz";
$som_vann = $o_folder.$t_id."_var_annotated.vcf.gz";
if(!$single_sample)	$som_vann = $o_folder.$t_id."-".$n_id."_var_annotated.vcf.gz";
$som_vqci = $o_folder.$t_id."_var_qci.vcf.gz";
if(!$single_sample)	$som_vqci = $o_folder.$t_id."-".$n_id."_var_qci.vcf.gz";
if (in_array("an", $steps))
{
	// annotate vcf into temp folder
	$tmp_folder1 = $parser->tempFolder();
	$tmp_vcf = $tmp_folder1."/".basename($som_gsvar, ".GSvar")."_var_annotated.vcf.gz";
	$parser->execTool("Pipelines/annotate.php", "-out_name ".basename($som_gsvar, ".GSvar")." -out_folder $tmp_folder1 -system $t_sys -vcf $som_v -t_col $t_id ".($single_sample?"":"-n_col $n_id")." -thres 8 -updown");
	if (!copy($tmp_vcf, $som_unfi))	trigger_error("Could not copy file '$tmp_vcf' to '$som_unfi'",E_USER_ERROR);

	// run somatic QC; this is run before project specific filtering since some qc parameters should not have additional filters set
	if(!$single_sample)
	{
		$t_bam = $t_folder.$t_id.".bam";
		$n_bam = $n_folder.$n_id.".bam";
		$links = array();
		if(is_file($t_folder.$t_id."_stats_fastq.qcML"))	$links[] = $t_folder.$t_id."_stats_fastq.qcML";
		if(is_file($t_folder.$t_id."_stats_map.qcML"))	$links[] = $t_folder.$t_id."_stats_map.qcML";
		if(is_file($n_folder.$n_id."_stats_fastq.qcML"))	$links[] = $n_folder.$n_id."_stats_fastq.qcML";
		if(is_file($n_folder.$n_id."_stats_map.qcML"))	$links[] = $n_folder.$n_id."_stats_map.qcML";
		$stafile3 = $o_folder.$t_id."-".$n_id."_stats_som.qcML";
		$parser->exec(get_path("ngs-bits")."SomaticQC","-tumor_bam $t_bam -normal_bam $n_bam ".(!empty($links)?"-links ".implode(" ",$links):"")." -somatic_vcf $som_unfi -target_bed ".$t_sys_ini['target_file']." -ref_fasta ".get_path("local_data")."/".$t_sys_ini['build'].".fa -out $stafile3",true);
	}

	// add project specific filters
	$extra = array();
	if($freebayes)	$extra[] = "-ignore_filter";
	if($t_sys_ini["target_file"]!="")	$extra[] = "-roi ".$t_sys_ini["target_file"];
	if(!$reduce_variants_filter)	$extra[] = "-keep";
	if($contamination > 0)	$extra[] = "-contamination $contamination";
	if($promoter!="")	$filter_set .= ",promoter -promoter $promoter";
	$parser->execTool("NGS/filter_vcf.php", "-in $som_unfi -out $som_vann -min_af $min_af -type $filter_set ".implode(" ", $extra));

	// annotate somatic NGSD-data
	// find processed sample with equal processing system for NGSD-annotation
	$extras = "-psname $t_id";
	if(is_valid_processingid($t_id)==false)
	{
		$extras = "";
		$tmp = get_processed_sample_name_by_processing_system($t_sys_ini['name_short'], false);
		if($tmp !== false)
		{
			$extras = "-psname $tmp";
			trigger_error("No valid processed sample id found ($t_id)! Used processing id $tmp instead. Statistics may be skewed!", E_USER_WARNING);
		}
	}

	// sort and dedup vcf comments
	$tmp = $parser->tempFile(".vcf");
	$s = Matrix::fromTSV($som_vann);
	$comments = $s->getComments();
	$details = get_processed_sample_info($t_id, false);
	$comments[] = gsvar_sample_header($t_id, array("IsTumor"=>"yes"), "#", "");
	if(!$single_sample)	$comments[] = gsvar_sample_header($n_id, array("IsTumor"=>"no"), "#", "");
	$s->setComments(sort_vcf_comments($comments));
	$s->toTSV($tmp);

	// annotate strand and family information for UID read groups
	if($t_sys_ini['type']=="Panel Haloplex HS")	$parser->exec(get_path("ngs-bits")."VariantAnnotateStrand", "-bam $t_bam -vcf $tmp -out $tmp  -hpHS ".substr($t_sys_ini['target_file'], 0, -4)."_amplicons.bed -name $t_id", true);
	if(!$single_sample && $n_sys_ini['type']=="Panel Haloplex HS")	$parser->exec(get_path("ngs-bits")."VariantAnnotateStrand", "-bam $n_bam -vcf $tmp -out $tmp  -hpHS ".substr($n_sys_ini['target_file'], 0, -4)."_amplicons.bed -name $n_id", true);
	if($t_sys_ini['type']=="Panel MIPs")	$parser->exec(get_path("ngs-bits")."VariantAnnotateStrand", "-bam $t_bam -vcf $tmp -out $tmp ".(isset($t_sys_ini['mip_file'])?$t_sys_ini['mip_file']:"/mnt/share/data/mipfiles/".$t_sys_ini["name_short"].".txt")." -name $t_id", true);
	if(!$single_sample && $n_sys_ini['type']=="Panel MIPs")	$parser->exec(get_path("ngs-bits")."VariantAnnotateStrand", "-bam $n_bam -vcf $tmp -out $tmp ".(isset($n_sys_ini['mip_file'])?$n_sys_ini['mip_file']:"/mnt/share/data/mipfiles/".$n_sys_ini["name_short"].".txt")." -name $n_id", true);

	// zip vcf file
	$parser->exec("bgzip", "-c $tmp > $som_vann", false);	// no output logging, because Toolbase::extractVersion() does not return
	$parser->exec("tabix", "-f -p vcf $som_vann", false);	// no output logging, because Toolbase::extractVersion() does not return

	// convert vcf to GSvar
	$extra = "-t_col $t_id ".($single_sample?"":"-n_col $n_id");
	if($t_sys_ini['type']=="Panel Haloplex HS" || starts_with($t_sys_ini['name_short'],"mi"))	$extra .= " -strand";
	$parser->execTool("NGS/vcf2gsvar_somatic.php", "-in $som_vann -out $som_gsvar $extra");

	// annotate NGSD and dbNFSP
	$parser->execTool("NGS/an_dbNFSPgene.php", "-in $som_gsvar -out $som_gsvar -build ".$t_sys_ini['build']);
	if (db_is_enabled("NGSD")) $parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode somatic",true);
	if (db_is_enabled("NGSD")) $parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode germline",true);
}

// qci
$som_vqci = $o_folder.$t_id."-".$n_id."_var_qci.vcf.gz";
if (in_array("ci", $steps))
{
	// add QCI output
	$parser->execTool("Tools/converter_vcf2qci.php", "-in $som_vann -out $som_vqci -pass");
}

// db
// add db import for qc parameters
// no import of variants (only db import for proper somatic pairs)
// import somatic variants to NGSD, check if sample within NGSD, otherwise skip QC import
if (in_array("db", $steps))
{
	if(is_valid_processingid($t_id) && (!$single_sample || is_valid_processingid($n_id)))
	{
		$parser->execTool("NGS/db_check_gender.php", "-in $t_bam -pid $t_id");
		if(!$single_sample)	$parser->execTool("NGS/db_check_gender.php", "-in $n_bam -pid $n_id");

		// import QC data tumor
		$log_db  = $t_folder."/".$t_id."_log4_db.log";
		$qcmls = $t_folder."/".$t_id."_stats_fastq.qcML ";
		if (is_file($t_folder."/".$t_id."_stats_map.qcML"))	$qcmls .= $t_folder."/".$t_id."_stats_map.qcML ";
		if (is_file($o_folder."/".$t_id."-".$n_id."_stats_som.qcML"))	$qcmls .= $o_folder.$t_id."-".$n_id."_stats_som.qcML ";		
		$parser->execTool("NGS/db_import_qc.php", "-id $t_id -files $qcmls -force -min_depth 0 --log $log_db");
		if (in_array("ma", $steps))	updateLastAnalysisDate($t_id, $t_bam);
		// check that tumor is flagged as tumor and normal not as tumor
		if(!isTumor($t_id))	trigger_error("Tumor $t_id is not flagged as tumor in NGSD",E_USER_WARNING);

		// import QC data normal
		if(!$single_sample)
		{
			$log_db  = $n_folder."/".$n_id."_log4_db.log";
			$qcmls  = $n_folder."/".$n_id."_stats_fastq.qcML ";
			if (is_file($n_folder."/".$n_id."_stats_map.qcML"))	$qcmls .= $n_folder."/".$n_id."_stats_map.qcML ";
			$parser->execTool("NGS/db_import_qc.php", "-id $n_id -files $qcmls -force -min_depth 0 --log $log_db");
			if (in_array("ma", $steps))	updateLastAnalysisDate($n_id, $t_bam);
			// check that tumor is flagged as tumor and normal not as tumor
			if(isTumor($n_id))	trigger_error("Normal $n_id is flagged as tumor in NGSD",E_USER_WARNING);
			// update normal entry for tumor
			if(updateNormalSample($t_id, $n_id))	trigger_error("Updated normal sample ($n_id) for tumor ($t_id) within NGSD.",E_USER_NOTICE);
		}

		// import variants (not for WGS) and if processing ids are valid
		if (!$single_sample && $t_sys_ini['type']!="WGS")	$parser->execTool("NGS/db_import_variants.php", "-id ".$t_id."-".$n_id." -var ".$som_gsvar." -build ".$t_sys_ini['build']." -force -mode somatic" /*--log $log_db"*/);
	}
	else	trigger_error("No DB import since no valid processing ID (T:".$t_id.($single_sample?"":"/N:".$n_id).")",E_USER_WARNING);
}