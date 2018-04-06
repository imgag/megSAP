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
$steps_all = array("ma", "vc", "an", "ci", "msi","db");
$parser->addString("steps", "Comma-separated list of processing steps to perform. Available are: ".implode(",", $steps_all), true, "ma,vc,an,msi,db");
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
$parser->addInfile("promoter","Bed file containing promoter regions. Will be used for filter column of vcf if filter 'not-promoter' is part of the filter_set.",true);
$parser->addEnum("clip", "Soft-clip overlapping read pairs.", true, array("sc","mfb","mfm","mfr","mfn"),"mfn");
$parser->addFlag("add_vc_folder", "Add folder containing variant calling results from variant caller.",true);
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));

// (0) preparations
// check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all))
	{
		trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	}
}

// tumor only or tumor normal
$single_sample = ($n_id === "na");
if ($single_sample)
{
	trigger_error("Single sample mode.", E_USER_NOTICE);
}
else
{
	trigger_error("Paired sample mode.", E_USER_NOTICE);
}

// check and generate output folders
$o_folder = $o_folder."/";
// tumor
$t_folder = $p_folder."/Sample_".$t_id."/";
$prefix_tum = $t_folder.$t_id;
$t_sys_ini = load_system($t_sys, $t_id);
if (!file_exists($t_folder))
{
	trigger_error("Tumor sample folder '$t_folder' does not exist.", E_USER_ERROR);
}

//amplicon mode
$amplicon = false;
if(!$t_sys_ini['shotgun'])
{
	$amplicon = true;
	trigger_error("Amplicon mode.",E_USER_NOTICE);
}

// normal
$n_folder = $p_folder."/Sample_".$n_id."/";
$prefix_nor = $n_folder.$n_id;

// get normal sample system and do some basic checks
if (!$single_sample)
{
	$n_sys_ini = load_system(($n_sys), $n_id);
	if (!file_exists($n_folder))
	{
		trigger_error("Reference sample folder '$n_folder' does not exist.", E_USER_ERROR);
	}

	// warn about missing target files
	if (empty($t_sys_ini['target_file']) || empty($n_sys_ini['target_file']))
	{
		trigger_error("Tumor system or normal system does not contain a target file.", E_USER_WARNING);
	}

	// warn when tumor and normal systems differ
	if ($t_sys_ini['name_short'] !== $n_sys_ini['name_short'])
	{
		trigger_error(sprintf("Tumor system '%s' and reference system '%s' are different!", $t_sys_ini['name_short'], $n_sys_ini['name_short']),
			E_USER_WARNING);
	}

	// fail in case of different genome builds
	if ($t_sys_ini['build'] !== $n_sys_ini['build'])
	{
		trigger_error(sprintf("Tumor system '%s' (%s) and reference system '%s' (%s) have different builds!",
			$t_sys_ini['name_short'], $t_sys_ini['build'], $n_sys_ini['name_short'], $n_sys_ini['build']),
			E_USER_ERROR);
	}
}


$target = $t_sys_ini['target_file'];

list($stdout) = $parser->exec(get_path("ngs-bits")."BedInfo", "-in $target", false);

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
		$args = [
			"-folder", $t_folder,
			"-name", $t_id,
			"-system", $t_sys,
			"--log", "{$t_folder}analyze_" . date('YmdHis',mktime()) . ".log",
			"-correction_n"
		];
		$parser->execTool("Pipelines/analyze.php", "-steps ma -no_abra ".implode(" ", $args));
	}

	// normal sample
	if (!$single_sample && !$smn)
	{
		$args = [
			"-folder", $n_folder,
			"-name", $n_id,
			"-system", $n_sys,
			"--log", "{$n_folder}analyze_" . date('YmdHis',mktime()) . ".log",
			"-correction_n"
		];
		$parser->execTool("Pipelines/analyze.php", "-steps ma,vc,an,cn -no_abra ".implode(" ", $args));
	}

	// overlap clipping, is not done by analyze.php to prevent import of QC data calculated on BamClipOverlap results
	if (!$no_softclip)
	{
		$clip_arg_map = [
			"mfb" => "-overlap_mismatch_baseq",
			"mfm" => "-overlap_mismatch_mapq",
			"mfr" => "-overlap_mismatch_remove",
			"mfn" => "-overlap_mismatch_basen",
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
	if ($abra && $t_sys_ini['type'] !== "WGS" && ($single_sample || $n_sys_ini['type'] !== "WGS"))
	{
		// intersect both target files
		$tmp_targets = $t_sys_ini['target_file'];
		if (!$single_sample)
		{
			$tmp_targets = $parser->tempFile("_intersected.bed");
			$commands_params = [];
			$commands_params[] = [get_path("ngs-bits")."BedIntersect", "-in ".$t_sys_ini['target_file']." -in2 ".$n_sys_ini['target_file']." -mode intersect"];
			$commands_params[] = [get_path("ngs-bits")."BedMerge", "-out $tmp_targets"];
			$parser->execPipeline($commands_params, "intersect target regions", true);
		}

		// indel realignment with ABRA
		$tmp1_t_bam = $parser->tempFile("_tumor.bam");
		$tmp1_n_bam = $parser->tempFile("_normal.bam");
		if (!$single_sample)
		{
			$tmp_in = "$n_bam $t_bam";
			$tmp_out = "$tmp1_n_bam $tmp1_t_bam";
		}
		else
		{
			$tmp_in = $t_bam;
			$tmp_out = $tmp1_t_bam;
		}

		$parser->execTool("NGS/indel_realign_abra.php", "-in $tmp_in -out $tmp_out -roi $tmp_targets -threads 8 -mer 0.02 -mad 5000");

		// copy realigned files to output folder and overwrite previous bam files
		$parser->moveFile($tmp1_t_bam, $t_bam);
		$parser->indexBam($t_bam, $threads);
		if (!$single_sample)
		{
			$parser->moveFile($tmp1_n_bam, $n_bam);
			$parser->indexBam($n_bam, $threads);
		}
	}
	
	// Low coverage statistics
	//get low_cov_statistics
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->exec(get_path("ngs-bits")."BedMerge", "-in $target -out $target_merged", true);
	//calculate low-coverage regions; 100x needed for 10 % sensitivity
	$td = 100;
	$low_cov = $o_folder."/".$t_id.($single_sample ? "" : "-".$n_id)."_stat_lowcov.bed";
	$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $target_merged -bam $t_bam -out $low_cov -cutoff ".$td, true);
	if(!$single_sample)
	{
		$nd = 100;
		$low_cov_n = $parser->tempFile("_nlowcov.bed");
		$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $target_merged -bam $n_bam -out $low_cov_n -cutoff ".$nd, true);
		$parser->exec(get_path("ngs-bits")."BedAdd", "-in $low_cov $low_cov_n -out $low_cov", true);
		$parser->exec(get_path("ngs-bits")."BedMerge", "-in $low_cov -out $low_cov", true);
	}
	//annotate gene names (with extended to annotate splicing regions as well)
	$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $low_cov -extend 25 -out $low_cov", true);
}

// check if samples are related
if (!$single_sample)
{
	if ($nsc)
	{
		trigger_error("Genotype correlation check has been disabled!", E_USER_WARNING);
	}
	else
	{
		$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in $t_bam $n_bam -mode bam -max_snps 4000", true);
		$correlation = explode("\t", $output[0][1])[3];
		if ($correlation < 0.8)
		{
			trigger_error("Genotype correlation of samples {$t_bam} and {$n_bam} is {$correlation}; should be at least 0.8!", E_USER_ERROR);
		}
	}
}

// variant calling
$prefix_som = $t_id . "-" . $n_id;		// "Tumor-Normal"
$prefix = $single_sample ? $o_folder . $t_id : $o_folder . $prefix_som;
$som_v     = $prefix . "_var.vcf.gz";	// variants
$som_cnv   = $prefix . "_cnvs.tsv";		// copy-number variants
$som_gsvar = $prefix . ".GSvar";		// GSvar variants
$som_bafs  = $prefix . "_bafs.igv";		// B-allele frequencies
$som_qci = $prefix . "_var_qci.vcf.gz";
$som_sv  = $prefix . "_var_structural.vcf.gz";	// structural variants (vcf)
$som_svt = $prefix . "_var_structural.tsv";		// structural variants (tsv)
$som_si  = $prefix . "_var_smallIndels.vcf.gz";	// small indels (manta)
if (in_array("vc", $steps))
{
	// structural variant calling
	// should be done before variant calling since strelka uses the smallIndel output
	if ($t_sys_ini['umi_type'] === "ThruPLEX")
	{
		trigger_error("Breakpoint detection deactivated for ThruPLEX samples.", E_USER_NOTICE);
	}
	else if ($amplicon)
	{
		trigger_error("Breakpoint detection deactivated for amplicon samples.", E_USER_NOTICE);
	}
	else
	{
		$args_manta = [
			"-t_bam", $t_bam,
			"-out", $som_sv,
			"-build", $t_sys_ini['build'],
			"-smallIndels", $som_si
		];

		if (!$single_sample)
		{
			$args_manta[] = "-bam $n_bam";
		}

		if ($t_sys_ini['type'] === "WES")
		{
			$args_manta[] = "-exome";
		}
		if ($add_vc_folder)
		{
			$args_manta[] = "-temp ".$o_folder."/manta";
		}
		if (!empty($t_sys_ini['target_file']))
		{
			$args_manta[] = "-target ".$t_sys_ini['target_file'];
		}
		$parser->execTool("NGS/vc_manta.php", implode(" ", $args_manta));
		$parser->execTool("Tools/converter_manta2tsv.php", "-in $som_sv -out $som_svt -tumor_id $t_id");
	}

	// variant calling
	if ($freebayes || $single_sample)	// variant calling using freebayes
	{
		// vc_freebayes uses ngs-bits that can not handle multi-sample vcfs
		$tmp1 = $parser->tempFile();
		$bams = $single_sample ? "$t_bam" : "$t_bam $n_bam";
		$parser->execTool("NGS/vc_freebayes.php","-bam $bams -out $tmp1 -build ".$t_sys_ini['build']." -min_af $min_af -target ".$t_sys_ini['target_file']);

		// find and rewrite header (tumor, normal sample columns)
		$tmp2 = $parser->tempFile();
		$s = Matrix::fromTSV($tmp1);
		$tmp_headers = $s->getHeaders();
		$tumor_col = NULL;
		$normal_col = NULL;
		for ($i=0; $i < count($tmp_headers); ++$i)
		{
			if (contains($tmp_headers[$i], $t_id))
			{
				$tumor_col = $i;
				$tmp_headers[$tumor_col] = $t_id;
			}
			else if (!$single_sample && contains($tmp_headers[$i],$n_id))
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
		for ($i=0; $i < $s->rows(); ++$i)
		{
			$row = $s->getRow($i);
			if ($row[$tumor_col]==".")
			{
				$s->set($i,$tumor_col,"./.:0:.:.:.:0:.,.,.");
			}
			if (!$single_sample && $row[$normal_col]==".")
			{
				$s->set($i,$normal_col,"./.:0:.:.:.:0:.,.,.");
			}
		}

		// add pedigree information for tumor and normal
		$pedigree = $single_sample ? "#PEDIGREE=<Tumor={$t_id}>" : "#PEDIGREE=<Tumor={$t_id},Normal={$n_id}>";
		$s->addComment($pedigree);

		// zip and index output file
		$s->toTSV($tmp2);
		$parser->exec("bgzip", "-c $tmp2 > $som_v", true);
		$parser->exec("tabix", "-f -p vcf $som_v", true);

		// filter for variants that are likely to be germline, included in strelka2
		$extra = array("-type somatic-lq","-keep");
		$parser->execTool("NGS/filter_vcf.php", "-in $som_v -out $som_v -min_af $min_af  ".implode(" ", $extra));
	}
	else
	{
		$args_strelka = [
			"-t_bam", $t_bam,
			"-n_bam", $n_bam,
			"-out", $som_v,
			"-build", $t_sys_ini['build'],
			"-target", $t_sys_ini['target_file']
		];
		if ($keep_all_variants_strelka)
		{
			$args_strelka[] = "-k";
		}
		if ($add_vc_folder)
		{
			$args_strelka[] = "-temp ".$o_folder."/strelka";
		}
		if (is_file($som_si))
		{
			$args_strelka[] = "-smallIndels $som_si";
		}
		$parser->execTool("NGS/vc_strelka2.php", implode(" ", $args_strelka));
	}

	// copy number variant calling
	$tmp_folder = $parser->tempFolder();

	// coverage for tumor sample
	$t_cov = "{$tmp_folder}/{$t_id}.cov";
	$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $t_bam -in ".$t_sys_ini['target_file']." -out $t_cov",true);

	// coverage for normal sample
	if (!$single_sample)
	{
		$n_cov = "{$tmp_folder}/{$n_id}.cov";
		$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $n_bam -in ".$n_sys_ini['target_file']." -out $n_cov", true);
	}

	// copy normal samplecoverage file to reference folder (has to be done before CnvHunter call to avoid analyzing the same sample twice)
	if (!$single_sample && is_valid_ref_sample_for_cnv_analysis($n_id))
	{
		//create reference folder if it does not exist
		$ref_folder = get_path("data_folder")."/coverage/".$n_sys_ini['name_short']."/";
		if (!is_dir($ref_folder))
		{
			mkdir($ref_folder);
		}
		//copy file
		$ref_file = $ref_folder.$n_id.".cov";
		$parser->copyFile($n_cov, $ref_file);
		$n_cov = $ref_file;
	}
	$ncov = $single_sample ? "" : "-n_cov $n_cov";
	$n = 20;
	if ($t_sys_ini['type'] == "WES")
	{
		$n = 30;
	}
	$parser->execTool("NGS/vc_cnvhunter.php", "-cov $t_cov $ncov -out $som_cnv -system $t_sys -min_corr 0 -seg $t_id -n $n");

	// add somatic baf file
	$baf_args = [
		"-in", $t_bam,
		"-out", $som_bafs,

	];
	if (!$single_sample)
	{
		$baf_args[] = "-n_in $n_bam";
	}
	if (!empty($t_sys_ini['target_file']))
	{
		$baf_args[] = "-target ".$t_sys_ini['target_file'];
	}
	$parser->execTool("NGS/mapping_baf.php", implode(" ", $baf_args));
}

// annotation
$som_unfi  = $prefix . "_var_annotated_unfiltered.vcf.gz";
$som_vann  = $prefix . "_var_annotated.vcf.gz";
$som_vqci  = $prefix . "_var_qci.vcf.gz";
if (in_array("an", $steps))
{
	// annotate vcf into temp folder
	$tmp_folder1 = $parser->tempFolder();
	$tmp_vcf = $tmp_folder1."/".$prefix_som."_var_annotated.vcf.gz";
	$parser->execTool("Pipelines/annotate.php",
		"-out_name $prefix_som -out_folder $tmp_folder1 -system $t_sys -vcf $som_v -t_col $t_id " . ($single_sample ? "" : "-n_col $n_id") . " -thres 8 -updown");
	$parser->copyFile($tmp_vcf, $som_unfi);

	// add donor annotation to full annotated vcf
	if (isset($donor_ids))
	{
		$donor_bams = [];
		foreach ($donor_ids as $donor_id)
		{
			$donor_bams[] = $p_folder."/Sample_{$donor_id}/{$donor_id}.bam";
		}
		$parser->execTool("NGS/vcf_somatic_donor.php", "-in_somatic {$som_unfi} -out_vcf {$som_unfi} -in_donor " . implode(" ", $donor_bams));
	}

	// run somatic QC; this is run before project specific filtering since some qc parameters should not have additional filters set
	if (!$single_sample)
	{
		$links = array_filter([
			$t_folder.$t_id."_stats_fastq.qcML",
			$t_folder.$t_id."_stats_map.qcML",
			$n_folder.$n_id."_stats_fastq.qcML",
			$n_folder.$n_id."_stats_map.qcML"
		], "file_exists");

		$stafile3 = $prefix . "_stats_som.qcML";
		$args_somaticqc = [
			"-tumor_bam", $t_bam,
			"-normal_bam", $n_bam,
			"-somatic_vcf", $som_unfi,
			"-target_bed", $t_sys_ini['target_file'],
			"-ref_fasta", get_path("local_data")."/".$t_sys_ini['build'].".fa",
			"-out", $stafile3
		];
		if (!empty($links))
		{
			$args_somaticqc[] = "-links";
			$args_somaticqc[] = implode(" ", $links);
		}

		$parser->exec(get_path("ngs-bits")."SomaticQC", implode(" ", $args_somaticqc), true);
	}

	// add project specific filters
	$filter_set = explode(",", $filter_set);
	if (isset($donor_ids))
	{
		$filter_set[] = "somatic-donor";
	}
	$args_filter_vcf = [
		"-in", $som_unfi,
		"-out", $som_vann,
		"-min_af", $min_af,
		"-type", implode(",", $filter_set)
	];
	// promoter BED file
	if (isset($promoter))
	{
		$args_filter_vcf[] = "-promoter $promoter";
	}
	if ($freebayes)
	{
		$args_filter_vcf[] = "-ignore_filter";
	}
	if ($t_sys_ini["target_file"] != "")
	{
		$args_filter_vcf[] = "-roi ".$t_sys_ini["target_file"];
	}
	if (!$reduce_variants_filter)
	{
		$args_filter_vcf[] = "-keep";
	}
	if ($contamination > 0)
	{
		$args_filter_vcf[] = "-contamination $contamination";
	}
	$parser->execTool("NGS/filter_vcf.php", implode(" ", $args_filter_vcf));

	// sort and dedup vcf comments
	$tmp = $parser->tempFile(".vcf");
	$s = Matrix::fromTSV($som_vann);
	$comments = $s->getComments();
	$details = get_processed_sample_info($t_id, false);
	$comments[] = gsvar_sample_header($t_id, array("IsTumor" => "yes"), "#", "");
	if (!$single_sample)
	{
		$comments[] = gsvar_sample_header($n_id, array("IsTumor" => "no"), "#", "");
	}
	$s->setComments(sort_vcf_comments($comments));
	$s->toTSV($tmp);

	// zip and index vcf file
	$parser->exec("bgzip", "-c $tmp > $som_vann", true);
	$parser->exec("tabix", "-f -p vcf $som_vann", true);

	// convert vcf to GSvar
	$extra = $single_sample ? "-t_col $t_id" : "-t_col $t_id -n_col $n_id";
	$parser->execTool("NGS/vcf2gsvar_somatic.php", "-in $som_vann -out $som_gsvar $extra");

	// annotate NGSD and dbNFSP
	$parser->execTool("NGS/an_dbNFSPgene.php", "-in $som_gsvar -out $som_gsvar -build ".$t_sys_ini['build']);
	if (db_is_enabled("NGSD"))
	{
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -psname $t_id -mode somatic", true);
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $som_gsvar -out $som_gsvar -psname $t_id -mode germline", true);
	}
}

// qci / CGI
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

			//icd10
			$laboratory_number = explode('_',$t_id)[0];
			$query = "SELECT ICD10DIAGNOSE,HPOTERM1 FROM `genlab8`.`v_ngs_sap` where labornummer = '$laboratory_number'";
			$result = $db->executeQuery($query);
			
			if(count($result) == 0) //"sometimes" GENLAB uses the full tumor ID as laboratory number
			{
				$query = "SELECT ICD10DIAGNOSE,HPOTERM1 FROM `genlab8`.`v_ngs_sap` where labornummer = '$t_id'";
				$result = $db->executeQuery($query);
			}
			
			if(count($result) == 0) //"sometimes" GENLAB uses the full tumor ID of the first processed sample as laboratory number
			{
				$laboratory_number = $laboratory_number."_01";
				$query = "SELECT ICD10DIAGNOSE,HPOTERM1 FROM `genlab8`.`v_ngs_sap` where labornummer = '$laboratory_number'";
				$result = $db->executeQuery($query);
			}
			
			$icd10_diagnosis = "";
			if(!empty($result))
			{
				$icd10_diagnosis = $result[0]['ICD10DIAGNOSE'];
			}
			
			
			//remove G and V from $icd10_diagnoses (sometimes there is a G or V assigned to its right)
			$icd10_diagnosis = rtrim($icd10_diagnosis,'G');
			$icd10_diagnosis = rtrim($icd10_diagnosis,'V');
			$hpo_diagnosis = "";
			if(!empty($result)) $hpo_diagnosis = $result[0]['HPOTERM1'];
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
						trigger_error("Could not assign CGI Cancer acronym uniquely to ICD10-diagnosis '$icd10_diagnosis'. Using standard cancertype CANCER instead.",E_USER_WARNING);
						$cancer_type = "CANCER";
					} else {
						$cancer_type = $hpo_cgi_dictionary[$hpo_diagnosis];
					}
				}
			} else {
				trigger_error("Could not assign CGI Cancer acronym to ICD10-diagnosis '$icd10_diagnosis'. Using standard cancertype CANCER instead.",E_USER_WARNING);
				$cancer_type = "CANCER";
			}
		}
	}

	$parameters = "";
	if(file_exists($som_vann))
	{
		$parameters = $parameters . " -mutations $som_vann";
	}
	if(file_exists($som_cnv))
	{
		$parameters = $parameters . " -cnas $som_cnv";
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

	//execution will not stop pipeline if it fails but display several errors
	$result_send_data = $parser->execTool("NGS/cgi_send_data.php", $parameters,false);
	$error_code = $result_send_data[2];
	if($error_code != 0)
	{
		trigger_error("step \"ci\" did not exit cleanly. Please check sample CGI files manually. cgi_send_data returned code ".$error_code,E_USER_WARNING);
	}
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
		$parameters = " -cnv_in $som_cnv -cnv_in_cgi $cgi_cnv_result_file -out $som_cnv";
		$parser->execTool("NGS/cgi_annotate_cnvs.php",$parameters);
	}
}

//msi: call microsatellite instabilities
if ($single_sample) trigger_error("Calling microsatellite instabilities is only possible for tumor normal pairs",E_USER_NOTICE);
else if (in_array("msi", $steps))
{
	$build = $n_sys_ini['build'];
	$reference_genome = get_path("data_folder") . "/genomes/" .$n_sys_ini['build'] . ".fa";
	
	//check whether file with loci exists in output folder
	//if not: intersect with loci file of reference
	$reference_loci_file = get_path("data_folder") . "/dbs/MANTIS/".$build."_msi_loci.bed";
	if(!file_exists($reference_loci_file))
	{
		trigger_error("Could not find loci reference file $reference_loci_file. Trying to generate it.",E_USER_WARNING);
		$parser->exec(get_path("mantis")."/tools/RepeatFinder","-i $reference_genome -o $reference_loci_file",true);
	}

	//file that contains MSI in target region -> is intersection of loci reference with target region
	$target_bed_file = $n_sys_ini['target_file'];
	$target_loci_file = $o_folder . "/" .$t_id. "-" . $n_id . "_msi_loci_target.bed";
	
	//target loci file is intersection of reference loci with target region
	if(!file_exists($target_loci_file));
	{
		$parameters = "-in ".$reference_loci_file." -in2 ".$target_bed_file ." -mode in -out ".$target_loci_file;
		$parser->exec(get_path("ngs-bits")."BedIntersect",$parameters,true);
	}
	$parameters = "-n_bam $n_bam -t_bam $t_bam -threads $threads -bed_file $target_loci_file -out $o_folder/mantis_test_output -build $reference_genome";
	$parser->execTool("NGS/detect_msi.php",$parameters);
}

// db
// add db import for qc parameters
// no import of variants (only db import for proper somatic pairs)
// import somatic variants to NGSD, check if sample within NGSD, otherwise skip QC import
if (in_array("db", $steps))
{
	if (!is_valid_processingid($t_id) || (!$single_sample && !is_valid_processingid($n_id)))
	{
		trigger_error("No database import since no valid processing ID (T:" . $t_id . ($single_sample ? "" : "/N:" . $n_id).")",
			E_USER_WARNING);
	}
	else
	{
		// check sex
		$parser->execTool("NGS/db_check_gender.php", "-in $t_bam -pid $t_id");

		// import qcML files
		$log_db = $prefix_tum."_log4_db.log";
		$qcmls = implode(" ", array_filter([
			$prefix_tum . "_stats_fastq.qcML",
			$prefix_tum . "_stats_map.qcML",
			$prefix . "_stats_som.qcML"
		], "file_exists"));
		$parser->execTool("NGS/db_import_qc.php", "-id $t_id -files $qcmls -force -min_depth 0 --log $log_db");

		// update last analysis date
		if (in_array("ma", $steps))
		{
			updateLastAnalysisDate($t_id, $t_bam);
		}

		// check tumor/normal flag
		if (!isTumor($t_id))
		{
			trigger_error("Tumor sample $t_id is not flagged as tumor in NGSD!", E_USER_WARNING);
		}

		// analogous steps for normal sample, plus additional
		if (!$single_sample)
		{
			// check sex
			$parser->execTool("NGS/db_check_gender.php", "-in $n_bam -pid $n_id");

			// import qcML files
			$log_db = $prefix_nor."_log4_db.log";
			$qcmls = implode(" ", array_filter([
				$prefix_nor . "_stats_fastq.qcML",
				$prefix_nor . "_stats_map.qcML",
				$prefix . "_stats_som.qcML"
			], "file_exists"));
			$parser->execTool("NGS/db_import_qc.php", "-id $n_id -files $qcmls -force -min_depth 0 --log $log_db");

			// update last analysis date
			if (in_array("ma", $steps))
			{
				updateLastAnalysisDate($n_id, $n_bam);
			}

			// check tumor/normal flag
			if (isTumor($n_id))
			{
				trigger_error("Normal sample $n_id is flagged as tumor in NGSD!", E_USER_WARNING);
			}

			// update normal sample entry for tumor, warn if different entry
			if (updateNormalSample($t_id, $n_id))
			{
				trigger_error("Updated normal sample ($n_id) for tumor ($t_id) in NGSD.", E_USER_NOTICE);
			}

			// import variants (not for WGS)
			if ($t_sys_ini['type'] !== "WGS")
			{
				$parser->execTool("NGS/db_import_variants.php", "-id {$t_id}-{$n_id} -var {$som_gsvar} -build ".$t_sys_ini['build']." -force -mode somatic");
			}
		}
	}
}

?>
