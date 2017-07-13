<?php

/**
	@page somatic_dna
		@TODO combine tumor normal pair and single sample mode
		@TODO recalculate QC-values at the end of the pipeline
		@TODO remove analyze.php and set up own combinations of tools
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
$parser->addString("filter_set","Filter set to use. Only if annotation step is selected. Multiple filters can be comma separated.",true,"non_coding_splicing,off_target,set_somatic");
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
$parser->addEnum("clip", "Soft-clip overlapping read pairs.", true, array("sc","mfb","mfm","mfr"),"sc");
$parser->addFlag("strelka1","Use strelka1 for variant calling.",true);
$parser->addFlag("add_vc_folder","Add folder containing variant calling results from variant caller.",true);
extract($parser->parse($argv));

// (0) preparations
// check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}
$tumor_only = ($n_id=="na");

// check and generate output folders
$o_folder = $o_folder."/";
$t_folder = $p_folder."/Sample_".$t_id."/";
if (!file_exists($t_folder))	trigger_error("Tumor-folder '$t_folder' does not exist.", E_USER_ERROR);
$n_folder = $p_folder."/Sample_".$n_id."/";
if (!$tumor_only && !file_exists($n_folder))	trigger_error("Reference-folder '$n_folder' does not exist.", E_USER_ERROR);
$t_sys_ini = load_system($t_sys, $t_id);

$amplicon = false;
if(!$t_sys_ini['shotgun'])
{
	$amplicon = true;
	trigger_error("Amplicon mode.",E_USER_NOTICE);
}

// (1) run analysis
$t_bam = $t_folder.$t_id.".bam";
$n_bam = $n_folder.$n_id.".bam";
// (1a) only tumor sample available	=>	freebayes
if($tumor_only)
{
	trigger_error("Single sample mode.",E_USER_NOTICE);

	$som_prefix = $t_id;
	$som_gsvar = $o_folder.$som_prefix.".GSvar";
	if (in_array("ma", $steps))
	{
		$args = "-steps ma";
		$parser->execTool("Pipelines/analyze.php", "-folder ".$t_folder." -out_folder ".$t_folder." -name $t_id -system $t_sys $args --log ".$t_folder."analyze_".date('YmdHis',mktime()).".log");

		if(!$no_softclip)
		{
			$extra = "";
			if($clip=="mfb")	$extra .= " -overlap_mismatch_baseq";
			if($clip=="mfm")	$extra .= " -overlap_mismatch_mapq";
			if($clip=="mfr")	$extra .= " -overlap_mismatch_remove";
//			if($amplicon)	$extra .= " -amplicon";
			
		
			$tmp1_t_bam = $parser->tempFile("_tumor.bam");
			$parser->exec(get_path("ngs-bits")."BamClipOverlap", " -in $t_bam -out $tmp1_t_bam $extra", true);
			$parser->exec(get_path("samtools"),"sort -T $tmp1_t_bam -o $t_bam $tmp1_t_bam", true);
			$parser->exec(get_path("ngs-bits")."BamIndex", "-in $t_bam", true);
		}
	}

	$som_v = $o_folder.$t_id."_var.vcf.gz";
	$som_cnv = $o_folder.$t_id."_var_cnvs.tsv";
	if (in_array("vc", $steps))
	{
		$tmp1 = $parser->tempFile();
		$parser->execTool("NGS/vc_freebayes.php", "-bam $t_bam -out $tmp1 -build ".$t_sys_ini['build']." -min_af ".$min_af." -target ".$t_sys_ini['target_file']);

		$tmp = $parser->tempFile();
		$s = Matrix::fromTSV($tmp1);
		// fix headers, replace id of sample column
		$tmp_headers = $s->getHeaders();
		for($i=0;$i<count($tmp_headers);++$i)
		{
			if(strpos($t_id,$tmp_headers[$i])===0)	$tmp_headers[$i] = $t_id;
		}
		$s->setHeaders($tmp_headers);
		// add tumor normal information, add pedigree information for SNPeff
		$s->addComment("#PEDIGREE=<Tumor=$t_id>");
		$s->toTSV($tmp);

		// zip and index output file
		$parser->exec("bgzip", "-c $tmp > $som_v", false);
		$parser->exec("tabix", "-f -p vcf $som_v", false);

		// copy number variant calling
		//copy coverage file to reference folder (has to be done before CnvHunter call to avoid analyzing the same sample twice)
		if (is_valid_ref_sample_for_cnv_analysis($n_id))
		{
			//create coverage file
			$tmp_folder = $parser->tempFolder();
			$cov_file = $tmp_folder."/{$n_id}.cov";
			$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $n_bam -in ".$n_sys_ini['target_file']." -out $cov_file", true);
			//create reference folder if it does not exist
			$ref_folder = get_path("data_folder")."/coverage/".$n_sys_ini['name_short']."/";
			if (!is_dir($ref_folder)) mkdir($ref_folder);
			
			//copy file
			$ref_file = $ref_folder.$n_id.".cov";
			copy2($cov_file, $ref_file);
			$cov_file = $ref_file;
		}
		$tmp_folder = $parser->tempFolder();
		$t_cov = $tmp_folder."/".basename($t_bam,".bam").".cov";
		$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $t_bam -in ".$t_sys_ini['target_file']." -out $t_cov",true);
		$parser->execTool("NGS/vc_cnvhunter.php", "-cov $t_cov -out $som_cnv -system $t_sys -min_reg 3 -min_corr 0 -seg $t_id");
	}

	$som_vcf = $o_folder.$t_id."_var_annotated.vcf.gz";
	$som_gsvar = $o_folder.$t_id.".GSvar";
	if (in_array("an", $steps))
	{
		// annotate vcf
		$tmp_folder1 = $parser->tempFolder();
		$parser->execTool("Pipelines/annotate.php", "-out_name ".basename($som_gsvar, ".GSvar")." -out_folder $tmp_folder1 -system $t_sys -vcf $som_v -t_col $t_id");

		// filter vcf to output folder
		$extra = array();
		if($t_sys_ini["target_file"]!="")	$extra[] = "-roi ".$t_sys_ini["target_file"];
		if(!$reduce_variants_filter)	$extra[] = "-keep";
		$parser->execTool("NGS/filter_vcf.php", "-in ${tmp_folder1}/".basename($som_gsvar, ".GSvar")."_var_annotated.vcf.gz -out $som_vcf -min_af $min_af -type $filter_set ".implode(" ", $extra));

		// sort vcf comments
		$tmp = $parser->tempFile(".vcf");
		$s = Matrix::fromTSV($som_vcf);
		$tmp_comments = sort_vcf_comments($s->getComments());
		$s->setComments($tmp_comments);
		$s->toTSV($tmp);
		
		// annotate strand if enrichment with UIDs is used
		if($t_sys_ini['type']=="Panel Haloplex HS")	$parser->exec(get_path("ngs-bits")."VariantAnnotateStrand", "-bam $t_bam -vcf $tmp -out $tmp  -hpHS ".substr($t_sys_ini['target_file'], 0, -4)."_amplicons.bed -name $t_id", true);
		if($t_sys_ini['type']=="Panel MIPs")
		{
			$mip_file = isset($t_sys_ini['mip_file']) ? $t_sys_ini['mip_file'] : "/mnt/share/data/mipfiles/".$t_sys_ini["name_short"].".txt";
			$parser->exec(get_path("ngs-bits")."VariantAnnotateStrand", "-bam $t_bam -vcf $tmp -out $tmp  -mip $mip_file -name $t_id", true);
		}
			
		// zip vcf file
		$parser->exec("bgzip", "-c $tmp > $som_vcf", false); //no output logging, because Toolbase::extractVersion() does not return
		$parser->exec("tabix", "-f -p vcf $som_vcf", false); //no output logging, because Toolbase::extractVersion() does not return	

		// convert vcf to GSvar
		$extra = "-t_col $t_id";
		if($t_sys_ini['type']=="Panel Haloplex HS" || $t_sys_ini['type']=="Panel MIPs")	$extra .= " -strand";
		$parser->execTool("NGS/vcf2gsvar_somatic.php", "-in $som_vcf -out $som_gsvar $extra");
		
		// dbNFSP
		$parser->execTool("NGS/an_dbNFSPgene.php", "-in $som_gsvar -out $som_gsvar");	
		
		// NGSD-data (somatic and germline)
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode somatic",true);
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode germline",true);
	}

	// add db import for qc parameters
	// no import of variants (only db import for proper somatic pairs)
	if (in_array("db", $steps) && is_valid_processingid($t_id))
	{
		$parser->execTool("NGS/db_check_gender.php", "-in $t_bam -pid $t_id");

		// import QC data tumor
		$log_db  = $t_folder."/".$t_id."_log4_db.log";
		$qc_fastq  = $t_folder."/".$t_id."_stats_fastq.qcML";
		$qc_map  = $t_folder."/".$t_id."_stats_map.qcML";
		$parser->execTool("NGS/db_import_qc.php", "-id $t_id -files $qc_fastq $qc_map -force -min_depth 0 --log $log_db");

		// update last_analysis date
		if(is_valid_processingid($t_id))
		{
			updateLastAnalysisDate($t_id, $t_bam);
		}
	}
}
// (1b) tumor normal pair => strelka or freebayes
else
{
	trigger_error("Tumor normal pair. Paired mode.",E_USER_NOTICE);

	// get ref system and do some basic checks
	$n_sys_ini = load_system(($n_sys), $n_id);
	if(empty($t_sys_ini['target_file']) || empty($n_sys_ini['target_file'])) trigger_error ("System tumor or system normal does not contain a target file.", E_USER_WARNING);
	if($t_sys_ini['name_short'] != $n_sys_ini['name_short']) trigger_error ("System tumor '".$t_sys_ini['name_short']."' and system reference '".$n_sys_ini['name_short']."' are different!", E_USER_WARNING);
	if($t_sys_ini['build'] != $n_sys_ini['build']) trigger_error ("System tumor '".$t_sys_ini['build']."' and system reference '".$n_sys_ini['build']."' do have different builds!", E_USER_ERROR);

	// map reference and tumor sample
	if (in_array("ma", $steps))
	{
		// mapping of tumor and reference sample, CAVE: no abra realignment and no soft-clipping, will be done later
		$args = "-steps ma -no_abra ";
		if(!$smt)	$parser->execTool("Pipelines/analyze.php", "-folder ".$t_folder." -name $t_id -system $t_sys ".$args." --log ".$t_folder."analyze_".date('YmdHis',mktime()).".log");
		if(!$smn)	$parser->execTool ("Pipelines/analyze.php", "-folder ".$n_folder." -name $n_id -system $n_sys ".$args." --log ".$n_folder."analyze_".date('YmdHis',mktime()).".log");

		if(!$no_softclip)
		{
			$extra = "";
			if($clip=="mfb")	$extra = "-overlap_mismatch_baseq";
			if($clip=="mfm")	$extra = "-overlap_mismatch_mapq";
			if($clip=="mfr")	$exta = "-overlap_mismatch_remove";
//			if($amplicon)	$extra .= " -amplicon";
		
			$tmp1_t_bam = $parser->tempFile("_tumor.bam");
			$parser->exec(get_path("ngs-bits")."BamClipOverlap", " -in $t_bam -out $tmp1_t_bam $extra", true);
			$parser->exec(get_path("samtools"),"sort -T $tmp1_t_bam -o $t_bam $tmp1_t_bam", true);
			$parser->exec(get_path("ngs-bits")."BamIndex", "-in $t_bam", true);
			
			$tmp1_n_bam = $parser->tempFile("_normal.bam");
			$parser->exec(get_path("ngs-bits")."BamClipOverlap", " -in $n_bam -out $tmp1_n_bam $extra", true);
			$parser->exec(get_path("samtools"),"sort -T $tmp1_t_bam -o $n_bam $tmp1_n_bam", true);
			$parser->exec(get_path("ngs-bits")."BamIndex", "-in $n_bam", true);
		}

		// indel realignment with ABRA
		if ($abra && $t_sys_ini['type']!="WGS" && $n_sys_ini['type']!="WGS")
		{
			// merge both target files
			$tmp1_targets = $parser->tempFile("_intersected.bed");
			$parser->exec(get_path("ngs-bits")."BedIntersect", "-in ".$t_sys_ini['target_file']." -in2 ".$n_sys_ini['target_file']." -out $tmp1_targets -mode intersect", true);
			$tmp2_targets = $parser->tempFile("_merged.bed");
			$parser->exec(get_path("ngs-bits")."BedMerge", "-in $tmp1_targets -out $tmp2_targets", true);

			// indel realignment with ABRA
			$tmp1_t_bam = $parser->tempFile("_tumor.bam");
			$tmp1_n_bam = $parser->tempFile("_normal.bam");
			$command = "php ".repository_basedir()."/src/NGS/indel_realign_abra.php -in $n_bam $t_bam -out $tmp1_n_bam $tmp1_t_bam -roi $tmp2_targets -threads 8 -mer 0.02 -mad 5000 2>&1";			
			$working_directory = realpath($p_folder);
			$parser->jobsSubmit(array($command), $working_directory, get_path("queues_high_mem"), true);
			
			// copy realigned files to output folder and overwrite previous bam files
			copy2($tmp1_n_bam, $n_bam);
			$parser->exec(get_path("ngs-bits")."BamIndex", "-in ".$n_bam, true);
			copy2($tmp1_t_bam, $t_bam);
			$parser->exec(get_path("ngs-bits")."BamIndex", "-in ".$t_bam, true);
		}
	}

	// check if samples are related
	if(!$nsc)
	{
		$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in1 $t_bam -in2 $n_bam -bam -max_snps 4000", true);
		$parts = explode(":", $output[0][1]);
		if ($parts[1]<0.8)	trigger_error("The genotype correlation of samples $t_bam and $n_bam is ".$parts[1]."; should be above 0.8!", E_USER_ERROR);
	}
		
	// variant calling
	$som_v = $o_folder.$t_id."-".$n_id."_var.vcf.gz";
	$som_qci = $o_folder.$t_id."-".$n_id."_var_qci.vcf.gz";
	$som_sv = $o_folder.$t_id."-".$n_id."_var_structural.vcf.gz";
	$som_si = $o_folder.$t_id."-".$n_id."_var_smallIndels.vcf.gz";
	$som_cnv = $o_folder.$t_id."-".$n_id."_cnvs.tsv";
	if (in_array("vc", $steps))
	{			
		// structural variant calling
		if($t_sys_ini['shotgun']==1)
		{
			$par = "";
			if($t_sys_ini['type']=="WES")	$par .= "-exome ";
			if($add_vc_folder)	$par .= "-temp ".dirname($som_v)."/variant_calling ";
			if(!$strelka1)	$par .= "-smallIndels $som_si ";
			$parser->execTool("NGS/vc_manta.php", "-t_bam $t_bam -bam $n_bam $par -out $som_sv -build ".$t_sys_ini['build']);
		}
		
		// variant calling
		if(!$freebayes)
		{
			$args = array();
			$args[] = "-build ".$t_sys_ini['build'];
			if ($keep_all_variants_strelka) $args[] = "-k";
			if ($amplicon) $args[] = "-amplicon";
			if($add_vc_folder)	$args[] = "-temp ".dirname($som_v)."/variant_calling";
			if(!$strelka1 && is_file($som_si))	$args[] = "-smallIndels $som_si";
			if(!$strelka1)	$parser->execTool("NGS/vc_strelka2.php", "-t_bam $t_bam -n_bam $n_bam -out $som_v ".implode(" ", $args));
			else 	$parser->execTool("NGS/vc_strelka.php", "-t_bam $t_bam -n_bam $n_bam -out $som_v ".implode(" ", $args));	//combined variant calling using strelka		
		}
		else	// combined variant calling using freebayes
		{
			// vc_freebayes uses owntools that can not handel multi-sample vcfs
			$tmp1 = $parser->tempFile();
			$parser->execTool("NGS/vc_freebayes.php", "-bam $t_bam $n_bam -out $tmp1 -build ".$t_sys_ini['build']." -min_af ".$min_af." -target ".$t_sys_ini['target_file']);

			// fix header
 			$tmp2 = $parser->tempFile();
			$s = Matrix::fromTSV($tmp1);
			$tmp_headers = $s->getHeaders();
			$tumor_col = NULL;
			$nomral_col = NULL;
			for($i=0;$i<count($tmp_headers);++$i)
			{
				if(strpos($t_id,$tmp_headers[$i])===0)
				{
					$tumor_col = $i;
					$tmp_headers[$tumor_col] = $t_id;
				}
				if(strpos($n_id,$tmp_headers[$i])===0)
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
				if($row[$normal_col]==".")
				{
					$s->set($i,$normal_col,"./.:0:.:.:.:0:.,.,.");
				}
			}
			// add tumor + normal information
			$s->addComment("#PEDIGREE=<Tumor=$t_id,Normal=$n_id>");	// add pedigree information for tumor and normal
			$s->toTSV($tmp2);
			// zip and index output file
			$parser->exec("bgzip", "-c $tmp2 > $som_v", false);	// no output logging, because Toolbase::extractVersion() does not return
			$parser->exec("tabix", "-f -p vcf $som_v", false);	// no output logging, because Toolbase::extractVersion() does not return
		}		
		
		// copy number variant calling
		$tmp_folder = $parser->tempFolder();
		$t_cov = $tmp_folder."/".basename($t_bam,".bam").".cov";
		$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $t_bam -in ".$t_sys_ini['target_file']." -out $t_cov",true);
		$n_cov = $tmp_folder."/".basename($n_bam,".bam").".cov";
		$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $n_bam -in ".$t_sys_ini['target_file']." -out $n_cov",true);
		$parser->execTool("NGS/vc_cnvhunter.php", "-cov $t_cov -n_cov $n_cov -out $som_cnv -system $t_sys -min_reg 3 -min_corr 0 -seg somatic");
	}

	// annotation
	$som_gsvar = $o_folder.$t_id."-".$n_id.".GSvar";
	$som_vann = $o_folder.$t_id."-".$n_id."_var_annotated.vcf.gz";
	$som_vqci = $o_folder.$t_id."-".$n_id."_var_qci.vcf.gz";
	if (in_array("an", $steps))
	{
		// annotate vcf into temp folder
		$tmp_folder1 = $parser->tempFolder();
		$parser->execTool("Pipelines/annotate.php", "-out_name ".basename($som_gsvar, ".GSvar")." -out_folder $tmp_folder1 -system $t_sys -vcf $som_v -t_col $t_id -n_col $n_id");

		// filter vcf to output folder
		$extra = array();
		if($t_sys_ini["target_file"]!="")	$extra[] = "-roi ".$t_sys_ini["target_file"];
		if(!$reduce_variants_filter)	$extra[] = "-keep";
		if($contamination > 0)	$extra[] = "-contamination $contamination";
		$parser->execTool("NGS/filter_vcf.php", "-in ${tmp_folder1}/".basename($som_gsvar, ".GSvar")."_var_annotated.vcf.gz -out $som_vann -min_af $min_af -type $filter_set ".implode(" ", $extra));

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
		$comments[] = gsvar_sample_header($n_id, array("IsTumor"=>"no"), "#", "");
		$s->setComments(sort_vcf_comments($comments));
		$s->toTSV($tmp);
			
		if($t_sys_ini['type']=="Panel Haloplex HS" || $t_sys_ini['type']=="Panel MIPs")
		{
			$parser->exec(get_path("ngs-bits")."VariantAnnotateStrand", "-bam $t_bam -vcf $tmp -out $tmp -name $t_id", true);
		}
		if($n_sys_ini['type']=="Panel Haloplex HS" || $n_sys_ini['type']=="Panel MIPs")
		{
			$parser->exec(get_path("ngs-bits")."VariantAnnotateStrand", "-bam $n_bam -vcf $tmp -out $tmp -name $n_id", true);
		}
		$parser->exec("bgzip", "-c $tmp > $som_vann", false);	// no output logging, because Toolbase::extractVersion() does not return
		$parser->exec("tabix", "-f -p vcf $som_vann", false);	// no output logging, because Toolbase::extractVersion() does not return

		// somatic QC
		$t_bam = $t_folder.$t_id.".bam";
		$n_bam = $n_folder.$n_id.".bam";
		$links = array($t_folder.$t_id."_stats_fastq.qcML",$t_folder.$t_id."_stats_map.qcML",$n_folder.$n_id."_stats_fastq.qcML",$n_folder.$n_id."_stats_map.qcML");
		$stafile3 = $o_folder.$t_id."-".$n_id."_stats_som.qcML";
		if(!$nsc)	$parser->exec(get_path("ngs-bits")."SomaticQC","-tumor_bam $t_bam -normal_bam $n_bam -links ".implode(" ",$links)." -somatic_vcf $som_vann -target_bed ".$t_sys_ini['target_file']." -out $stafile3",true);

		// convert vcf to GSvar
		$extra = "-t_col $t_id -n_col $n_id";
		if($t_sys_ini['type']=="Panel Haloplex HS" || strpos($t_sys_ini['name_short'],"mi")===0)	$extra .= " -strand";
		$parser->execTool("NGS/vcf2gsvar_somatic.php", "-in $som_vann -out $som_gsvar $extra");
		$parser->execTool("NGS/an_dbNFSPgene.php", "-in $som_gsvar -out $som_gsvar -build ".$t_sys_ini['build']);
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode somatic",true);
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode germline",true);
	}

	// qci
	$som_vqci = $o_folder.$t_id."-".$n_id."_var_qci.vcf.gz";
	if (in_array("ci", $steps))
	{
		// add QCI output
		$parser->execTool("Tools/converter_vcf2qci.php", "-in $som_vann -out $som_vqci -pass");
	}

	// import somatic variants to NGSD, check if sample within NGSD, otherwise skip QC import
	if (in_array("db", $steps))
	{
		if(is_valid_processingid($t_id) && is_valid_processingid($n_id))
		{
			$parser->execTool("NGS/db_check_gender.php", "-in $t_bam -pid $t_id");
			$parser->execTool("NGS/db_check_gender.php", "-in $n_bam -pid $n_id");

			// import QC data tumor
			$log_db  = $t_folder."/".$t_id."_log4_db.log";
			$qcmls = $t_folder."/".$t_id."_stats_fastq.qcML ";
			if (in_array("ma", $steps))	$qcmls .= $t_folder."/".$t_id."_stats_map.qcML ";
			if (!$nsc && !$tumor_only && in_array("an", $steps))	$qc_som  = $o_folder.$t_id."-".$n_id."_stats_som.qcML ";
			$parser->execTool("NGS/db_import_qc.php", "-id $t_id -files $qcmls -force -min_depth 0 --log $log_db");

			// import QC data normal
			$log_db  = $n_folder."/".$n_id."_log4_db.log";
			$qcmls  = $n_folder."/".$n_id."_stats_fastq.qcML ";
			if (in_array("ma", $steps))	$qcmls .= $n_folder."/".$n_id."_stats_map.qcML ";
			$parser->execTool("NGS/db_import_qc.php", "-id $n_id -files $qcmls -force -min_depth 0 --log $log_db");

			// update last_analysis date
			if (in_array("ma", $steps))	updateLastAnalysisDate($t_id, $t_bam);
			if (in_array("ma", $steps))	updateLastAnalysisDate($n_id, $t_bam);
			
			// check that tumor is flagged as tumor and normal not as tumor
			if(!isTumor($t_id))	trigger_error("Tumor $t_id is not flagged as tumor in NGSD",E_USER_WARNING);
			if(isTumor($n_id))	trigger_error("Normal $n_id is flagged as tumor in NGSD",E_USER_WARNING);
			
			// update normal entry for tumor
			if(updateNormalSample($t_id, $n_id))	trigger_error("Updated normal sample ($n_id) for tumor ($t_id) within NGSD.",E_USER_NOTICE);
			
			// import variants (not for WGS) and if processing ids are valid
			if ($t_sys_ini['type']!="WGS" && $n_sys_ini['type']!="WGS")	$parser->execTool("NGS/db_import_variants.php", "-id ".$t_id."-".$n_id." -var ".$som_gsvar." -build ".$t_sys_ini['build']." -force -mode somatic" /*--log $log_db"*/);
		}
		else	trigger_error("No DB import since no valid processing ID (T:".$t_id."/N:".$n_id.")",E_USER_WARNING);
	}
}