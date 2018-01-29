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
$parser->addString("cancer_type","Tumor type. See CancerGenomeInterpreter.org for nomenclature. If not set, megSAP will try to resolve cancer type from GENLAB.",true);
$parser->addString("filter_set","Filter set to use. Only if annotation step is selected. Multiple filters can be comma separated.",true,"synonymous,not-coding-splicing");
// optional
$parser->addStringArray("donor_ids", "Donor sample IDs in case of transplanted patients.", true);
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
	$no_softclip = true;
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

		$parser->execTool("NGS/indel_realign_abra.php", "-in $tmp_in -out $tmp_out -roi $tmp_targets -threads 8 -mer 0.02 -mad 5000");
		
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
$som_bafs = $o_folder.$t_id."-".$n_id."_bafs.igv";
if (in_array("vc", $steps))
{
	// structural variant calling, should be done before variant calling since strelka uses the smallIndel output
	if(!$single_sample)
	{
		if($t_sys_ini['umi_type']!="ThruPLEX")
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
		else trigger_error("Breakpoint detection deactivated for ThruPLEX samples.",E_USER_NOTICE);
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
	$args = array();
	if(!$single_sample)	$args[] = "-n_in $n_bam";
	if(!empty($t_sys['target_file']))	$args[] = "-target ".$t_sys['target_file'];
	$parser->execTool("NGS/mapping_baf.php", "-in $t_bam -out $som_bafs ".implode(" ",$args));
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
	$parser->copyFile($tmp_vcf, $som_unfi);
	
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

	// add donor annotation to full annotated vcf
	if (isset($donor_ids))
	{
		$donor_bams = [];
		foreach($donor_ids as $donor_id)
		{
			$donor_bams[] = $p_folder."/Sample_{$donor_id}/{$donor_id}.bam";
		}
		$parser->execTool("NGS/vcf_somatic_donor.php", "-in {$som_unfi} -out {$som_unfi} -in_donor " . implode(" ", $donor_bams));
	}

	// add project specific filters
	$extra = array();
	if($freebayes)	$extra[] = "-ignore_filter";
	if($t_sys_ini["target_file"]!="")	$extra[] = "-roi ".$t_sys_ini["target_file"];
	if(!$reduce_variants_filter)	$extra[] = "-keep";
	if($contamination > 0)	$extra[] = "-contamination $contamination";
	// add somatic-donor filter if donor samples are specified
	if (isset($donor_ids))	$filter_set .= ",somatic-donor";
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

	// zip vcf file
	$parser->exec("bgzip", "-c $tmp > $som_vann", false);	// no output logging, because Toolbase::extractVersion() does not return
	$parser->exec("tabix", "-f -p vcf $som_vann", false);	// no output logging, because Toolbase::extractVersion() does not return

	// convert vcf to GSvar
	$extra = "-t_col $t_id ".($single_sample?"":"-n_col $n_id");
	$parser->execTool("NGS/vcf2gsvar_somatic.php", "-in $som_vann -out $som_gsvar $extra");

	// annotate NGSD and dbNFSP
	$parser->execTool("NGS/an_dbNFSPgene.php", "-in $som_gsvar -out $som_gsvar -build ".$t_sys_ini['build']);
	if (db_is_enabled("NGSD")) $parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode somatic",true);
	if (db_is_enabled("NGSD")) $parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -mode germline",true);
}

// qci / CGI
//file names, file path is expected in o_folder
$som_vqci = $o_folder.$t_id."_var_qci.vcf.gz";
if(!$single_sample) $som_vqci = $o_folder.$t_id."-".$n_id."_var_qci.vcf.gz";
$som_vann = $o_folder.$t_id."_var_annotated.vcf.gz";
if(!$single_sample)	$som_vann = $o_folder.$t_id."-".$n_id."_var_annotated.vcf.gz";
$som_cnvs = $o_folder.$t_id."_cnvs.tsv";
if(!$single_sample) $som_cnvs = $o_folder.$t_id."-".$n_id."_cnvs.tsv";
$som_gsvar = $o_folder.$t_id.".GSvar";
if(!$single_sample)	$som_gsvar = $o_folder.$t_id."-".$n_id.".GSvar";
//TODO: implementation for translocation files
if (in_array("ci", $steps))
{
	// add QCI output
	$parser->execTool("Tools/converter_vcf2qci.php", "-in $som_vann -out $som_vqci -pass");

	//check whether genlab credentials are available
	$db_hosts = get_path("db_host",false);
	if(!array_key_exists("GL8",$db_hosts))
	{
		trigger_error("No credentials for Genlab db found. Using generic cancertype CANCER",E_USER_WARNING);
		$cancer_type = "CANCER";
	}
	else
	{
		//if no cancer_type is set try to resolve cancer_type from GENLAB
		if(!isset($cancer_type))
		{
			//database connection to GENLAB
			$db = DB::getInstance("GL8");
			$laboratory_number = explode('_',$t_id)[0];
			$query = "SELECT ICD10DIAGNOSE,HPOTERM1 FROM `genlab8`.`v_ngs_sap` where labornummer = '$laboratory_number'";
			$result = $db->executeQuery($query);

			//icd10
			$icd10_diagnosis = "";
			if(!empty($result)) $icd10_diagnosis = $result[0]['ICD10DIAGNOSE'];
			//remove G and V from $icd10_diagnoses (sometimes there is a G or V assigned to its right)
			$icd10_diagnosis = rtrim($icd10_diagnosis,'G');
			$icd10_diagnosis = rtrim($icd10_diagnosis,'V');
			
			$hpo_diagnosis = $result[0]['HPOTERM1'];
			if(empty($icd10_diagnosis))
			{
				trigger_error("There is no ICD10 diagnosis set in Genlab.",E_USER_WARNING);
			}
			
			$dictionary = Matrix::fromTSV(repository_basedir()."/data/dbs/Ontologies/icd10_cgi_dictionary.tsv");
			$words_diagnoses = $dictionary->getCol($dictionary->getColumnIndex("icd10_code"));
			$words_cgi_acronyms = $dictionary->getCol($dictionary->getColumnIndex("cgi_acronym"));
			$words_hpo_terms = $dictionary->getCol($dictionary->getColumnIndex("hpo_term"));
			$icd10_cgi_dictionary = array_combine($words_diagnoses,$words_cgi_acronyms);
			$cgi_hpo_dictionary = array_combine($words_diagnoses,$words_hpo_terms);
			//if there is a cancer acronym annotated to ICD10 diagnosis set cancertype according to ICD 10 diagnosis
			if(!empty($icd10_cgi_dictionary[$icd10_diagnosis]))
			{
				//check whether assignment of ICD10 to CGI acronyms is unique
				$cgi_acronyms = explode(',',$icd10_cgi_dictionary[$icd10_diagnosis]);
				$hpo_terms = explode(',',$cgi_hpo_dictionary[$icd10_diagnosis]);
				if(count($cgi_acronyms) == 1) //unique assignment
				{
					$cancer_type = $icd10_cgi_dictionary[$icd10_diagnosis];
				} 
				else //not unique: try to assign a CGI cancer acronym via HPO terms
				{
					$hpo_cgi_dictionary = array_combine($hpo_terms,$cgi_acronyms);

					if(!array_key_exists($hpo_diagnosis,$hpo_cgi_dictionary))
					{
						trigger_error("Could not assign CGI Cancer acronym uniquely to ICD10-diagnosis $icd10_diagnosis. Using standard cancertype CANCER instead.",E_USER_WARNING);
						$cancer_type = "CANCER";
					} else {
						$cancer_type = $hpo_cgi_dictionary[$hpo_diagnosis];
					}
				}
			} else {
				trigger_error("Could not assign CGI Cancer acronym to ICD10-diagnosis $icd10_diagnosis. Using standard cancertype CANCER instead.",E_USER_WARNING);
				$cancer_type = "CANCER";
			}
		}
	}

	$parameters = "";
	if(file_exists($som_vann))
	{
		$parameters = $parameters . " -mutations $som_vann";
	}
	if(file_exists($som_cnvs))
	{
		$parameters = $parameters . " -cnas $som_cnvs";
	}
	if($cancer_type != "")
	{
		$parameters = $parameters . " -cancertype $cancer_type";
	}
	$parameters = $parameters . " -o_folder $o_folder";

	//if we know genes in target region we set parameter for this file
	$genes_in_target_region =  dirname($t_sys_ini["target_file"]) . "/" .basename($t_sys_ini["target_file"],".bed")."_genes.txt";
	if(file_exists($genes_in_target_region))
	{
		$parameters = $parameters . " -t_region $genes_in_target_region";
	}
	
	$parser->execTool("NGS/cgi_send_data.php", $parameters);
	
	$parameters = "";
	
	//only try annotate SNVS to GSVar file if $som_vann (variant file) was uploaded to CGI
	$cgi_snv_result_file = $o_folder . "/" .$t_id ."-".$n_id . "_cgi_mutation_analysis.tsv";
	if(file_exists($cgi_snv_result_file) && file_exists($som_gsvar))
	{
		$parameters = " -gsvar_in $som_gsvar -cgi_snv_in $cgi_snv_result_file -out $som_gsvar";
		$parser->execTool("NGS/cgi_snvs_to_gsvar.php",$parameters);	
	}
	//annotate CGI cnv genes to cnv input file (which was originally created by CNVHunter)
	$cgi_cnv_result_file = $o_folder . "/" .$t_id ."-".$n_id . "_cgi_cnv_analysis.tsv";
	if(file_exists($cgi_cnv_result_file))
	{
		$parameters = " -cnv_in $som_cnvs -cnv_in_cgi $cgi_cnv_result_file -out $som_cnvs";
		$parser->execTool("NGS/cgi_annotate_cnvs.php",$parameters);
	}
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