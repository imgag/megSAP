<?php

/**
 * @page somatic_dna
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_dna", "Differential analysis of tumor and normal DNA samples.");

//mandatory
$parser->addInfile("t_bam", "Tumor sample BAM file.", false);
$parser->addString("out_folder", "Output folder.", false);

//optional
$parser->addInfile("n_bam", "Normal sample BAM file.", true);
$parser->addInfile("t_rna_bam", "Tumor RNA sample BAM file.", true);

$parser->addString("prefix", "Output file prefix.", true, "somatic");

$steps_all = array("vc", "cn", "an", "ci", "msi","db");
$parser->addString("steps", "Comma-separated list of steps to perform:\n" .
	"vc=variant calling, an=annotation, ci=CGI annotation,\n" .
	"cn=copy-number analysis, msi=microsatellite analysis, db=database import",
	true, "vc,cn,an,ci,msi,db");

$parser->addString("cancer_type", "Tumor type, see CancerGenomeInterpreter.org for nomenclature (resolved from GENLAB if not set).", true);

$parser->addInfile("system",  "Processing system file used for tumor DNA sample (resolved from NGSD via tumor BAM by default).", true);
$parser->addInfile("n_system",  "Processing system file used for normal DNA sample (resolved from NGSD via tumor BAM by default).", true);

$parser->addFlag("skip_correlation", "Skip sample correlation check.");

//TODO remove this
$parser->addFlag("include_germline", "Include germline variant annotation.");

//default cut-offs
$parser->addFloat("min_af", "Allele frequency detection limit.", true, 0.05);
$parser->addFloat("min_correlation", "Minimum correlation for tumor/control pair.", true, 0.8);
$parser->addFloat("min_depth_t", "Tumor sample coverage cut-off for low coverage statistic.", true, 100);
$parser->addFloat("min_depth_n", "Control sample coverage cut-off for low coverage statistic.", true, 100);
$parser->addInt("min_cov_files", "Minimum number of required coverage files for CNV calling.", true, 20);

$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all))
	{
		trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	}
}

//tumor only or tumor-normal
$single_sample = !isset($n_bam);

//output prefix
$full_prefix = "{$out_folder}/{$prefix}";

//IDs, system and target region
$t_id = basename($t_bam, ".bam");
$n_id = basename($n_bam, ".bam");
$sys = load_system($system, $t_id);
if (!$single_sample)
{
	$n_sys = load_system($n_system, $n_id);
};
$roi = $sys["target_file"];

//sample similiarty check
$bams = array_filter([$t_bam, $n_bam, $t_rna_bam]);
if (count($bams) > 1)
{
    if ($skip_correlation)
    {
        trigger_error("Genotype correlation check has been disabled!", E_USER_WARNING);
    }
    else
    {
    	$args_similarity = [
    		"-in", implode(" ", $bams),
			"-mode", "bam",
			"-max_snps", 4000
		];
		if (!empty($roi))
		{
			$args_similarity[] = "-roi";
			$args_similarity[] = $roi;
		}
        $output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", implode(" ", $args_similarity), true);

		//extract columen 3 from output
		$table = array_map(
			function($str) { return explode("\t", $str); },
			array_slice($output[0], 1)
		);
		$correlation = array_column($table, 3);
        if (min($correlation) < $min_correlation)
        {
            trigger_error("Genotype correlation lower than {$min_correlation}!\n" . implode("\n", $output[0]), E_USER_ERROR);
        }
    }
}

//low coverage statistics
$low_cov = "{$full_prefix}_stat_lowcov.bed";					// low coverage BED file
if ($sys['type'] !== "WGS" && !empty($roi))
{
	$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $roi -bam $t_bam -out $low_cov -cutoff $min_depth_t", true);
	//combined tumor and normal low coverage files
	//normal coverage is calculated only for tumor target region
	if(!$single_sample)
	{
		$low_cov_n = $parser->tempFile("_nlowcov.bed");
		$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $roi -bam $n_bam -out $low_cov_n -cutoff $min_depth_n", true);
		$parser->execPipeline([
			[get_path("ngs-bits")."BedAdd", "-in $low_cov $low_cov_n"],
			[get_path("ngs-bits")."BedMerge", "-out $low_cov"]
		], "merge low coverage BED files");
	}
	// annotate with gene names
	if (db_is_enabled("NGSD"))
	{
		$parser->exec(get_path("ngs-bits") . "BedAnnotateGenes", "-in $low_cov -extend 25 -out $low_cov", true);
	}
}

//variant calling
$manta_indels = $full_prefix . "_var_smallIndels.vcf.gz";		// small indels from manta
$manta_sv     = $full_prefix . "_var_structural.vcf.gz";		// structural variants (vcf)
$manta_sv_tsv = $full_prefix . "_var_structural.tsv";			// structural variants (tsv)
$variants     = $full_prefix . "_var.vcf.gz";					// variants
$ballele      = $full_prefix . "_bafs.igv";						// B-allele frequencies
if (in_array("vc", $steps))
{
	// structural variant calling
	if (!$sys['shotgun'])
	{
		trigger_error("Structural variant calling deactivated for amplicon samples.", E_USER_NOTICE);
	}
	else if ($sys['umi_type'] === "ThruPLEX")
	{
		trigger_error("Structural variant calling deactivated for ThruPLEX samples.", E_USER_NOTICE);
	}
	else
	{
		$args_manta = [
			"-t_bam", $t_bam,
			"-out", $manta_sv,
			"-build", $sys['build'],
			"-smallIndels", $manta_indels,
			"-threads", $threads
		];
		if (!$single_sample)
		{
			$args_manta[] = "-bam $n_bam";
		}
		if ($sys['type'] === "WES")
		{
			$args_manta[] = "-exome";
		}
		if (!empty($roi))
		{
			$args_manta[] = "-target {$roi}";
		}
		$parser->execTool("NGS/vc_manta.php", implode(" ", $args_manta));
		$parser->execTool("Tools/converter_manta2tsv.php", "-in $manta_sv -out $manta_sv_tsv -tumor_id $t_id");
	}

	// variant calling
	if ($single_sample)	// variant calling using freebayes
	{
		// vc_freebayes uses ngs-bits that can not handle multi-sample vcfs
		$tmp1 = $parser->tempFile();
		$bams = $single_sample ? "$t_bam" : "$t_bam $n_bam";
		$parser->execTool("NGS/vc_freebayes.php", "-bam $bams -out $tmp1 -build ".$sys['build']." -min_af $min_af -target ".$roi);

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

		// zip and index output file
		$s->toTSV($tmp2);
		$parser->exec("bgzip", "-c $tmp2 > $variants", true);
		$parser->exec("tabix", "-f -p vcf $variants", true);
	}
	else
	{
		$args_strelka = [
			"-t_bam", $t_bam,
			"-n_bam", $n_bam,
			"-out", $variants,
			"-build", $sys['build']
		];
		if (!empty($roi))
		{
			$args_strelka[] = "-target $roi";
		}
		if (is_file($manta_indels))
		{
			$args_strelka[] = "-smallIndels $manta_indels";
		}
		$parser->execTool("NGS/vc_strelka2.php", implode(" ", $args_strelka));
	}

	// add somatic BAF file
	if (in_array($sys['build'], [ "hg19", "GRCh37" ]))
	{
	$baf_args = [
		"-in {$t_bam}",
		"-out {$ballele}", 
		"-build ".$sys['build']
	];
	if (!$single_sample)
	{
		$baf_args[] = "-n_in $n_bam";
	}
	if (!empty($roi))
	{
		$baf_args[] = "-target {$roi}";
	}
	$parser->execTool("NGS/mapping_baf.php", implode(" ", $baf_args));
}

//CNV calling
$cnvs = $full_prefix . "_cnvs.tsv";								// copy-number variants
if (in_array("cn", $steps) && $sys['type'] !== "WGS" && !empty($sys['target_file']) && ($single_sample || !empty($n_sys['target_file'])))
{
	// copy number variant calling
	$tmp_folder = $parser->tempFolder();

	// coverage for tumor sample
	$t_cov = "{$tmp_folder}/{$t_id}.cov";
	$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $t_bam -in $roi -out $t_cov",true);

	// coverage for normal sample
	if (!$single_sample)
	{
		$n_cov = "{$tmp_folder}/{$n_id}.cov";
		$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $n_bam -in ".$n_sys['target_file']." -out $n_cov", true);
	}

	// copy normal samplecoverage file to reference folder (has to be done before CnvHunter call to avoid analyzing the same sample twice)
	if (!$single_sample && is_valid_ref_sample_for_cnv_analysis($n_id))
	{
		//create reference folder if it does not exist
		$ref_folder = get_path("data_folder")."/coverage/".$n_sys['name_short']."/";
		create_directory($ref_folder);
		//copy file
		$ref_file = $ref_folder.$n_id.".cov";
		$parser->copyFile($n_cov, $ref_file);
		$n_cov = $ref_file;
	}
	$ncov = $single_sample ? "" : "-n_cov $n_cov";
	//at least 30 coverage files for WES
	if ($sys['type'] == "WES")
	{
		$min_cov_files = max($min_cov_files, 30);
	}
	$parser->execTool("NGS/vc_cnvhunter.php", "-cov $t_cov $ncov -out $cnvs -system $system -min_corr 0 -seg $t_id -n $min_cov_files");
}

//annotation
$variants_annotated = $full_prefix . "_var_annotated.vcf.gz";	// annotated variants
$variants_gsvar     = $full_prefix . ".GSvar";					// GSvar variants
$somaticqc          = $full_prefix . "_stats_som.qcML";			// SomaticQC qcML
if (in_array("an", $steps))
{
	// annotate vcf into temp folder
	$tmp_folder1 = $parser->tempFolder();
	$tmp_vcf = "{$tmp_folder1}/{$prefix}_var_annotated.vcf.gz";
	$parser->execTool("Pipelines/annotate.php", "-out_name $prefix -out_folder $tmp_folder1 -system $system -vcf $variants -somatic -updown");

	// add donor annotation to full annotated vcf
	if (isset($donor_ids))
	{
		$donor_bams = [];
		foreach ($donor_ids as $donor_id)
		{
			//TODO fix
			$donor_bams[] = $p_folder."/Sample_{$donor_id}/{$donor_id}.bam";
		}
		$parser->execTool("NGS/vcf_somatic_donor.php", "-in_somatic {$tmp_vcf} -out_vcf {$tmp_vcf} -in_donor " . implode(" ", $donor_bams));
	}

	// run somatic QC
	if (!$single_sample)
	{
		$links = array_filter([
			dirname($t_bam)."/{$t_id}_stats_fastq.qcML",
			dirname($t_bam)."/{$t_id}_stats_map.qcML",
			dirname($n_bam)."/{$n_id}_stats_fastq.qcML",
			dirname($n_bam)."/{$n_id}_stats_map.qcML"
		], "file_exists");

		$args_somaticqc = [
			"-tumor_bam", $t_bam,
			"-normal_bam", $n_bam,
			"-somatic_vcf", $tmp_vcf,
			"-target_bed", $roi,
			"-ref_fasta", get_path("local_data")."/".$sys['build'].".fa",
			"-out", $somaticqc
		];
		if (!empty($links))
		{
			$args_somaticqc[] = "-links";
			$args_somaticqc[] = implode(" ", $links);
		}

		$parser->exec(get_path("ngs-bits")."SomaticQC", implode(" ", $args_somaticqc), true);
	}

	// sort and dedup vcf comments
	$tmp = $parser->tempFile(".vcf");
	$s = Matrix::fromTSV($tmp_vcf);
	$comments = $s->getComments();
	$comments[] = gsvar_sample_header($t_id, array("IsTumor" => "yes"), "#", "");
	if (!$single_sample)
	{
		$comments[] = gsvar_sample_header($n_id, array("IsTumor" => "no"), "#", "");
	}
	$s->setComments(sort_vcf_comments($comments));
	$s->toTSV($tmp);

	// zip and index vcf file
	$parser->exec("bgzip", "-c $tmp > $variants_annotated", true);
	$parser->exec("tabix", "-f -p vcf $variants_annotated", true);

	// convert vcf to GSvar
	$extra = $single_sample ? "-t_col $t_id" : "-t_col $t_id -n_col $n_id";
	$parser->execTool("NGS/vcf2gsvar_somatic.php", "-in $variants_annotated -out $variants_gsvar $extra");

	// annotate NGSD and dbNFSP
	$parser->execTool("NGS/an_dbNFSPgene.php", "-in $variants_gsvar -out $variants_gsvar -build ".$sys['build']);
	if (db_is_enabled("NGSD"))
	{
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $variants_gsvar -out $variants_gsvar -psname $t_id -mode somatic", true);
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $variants_gsvar -out $variants_gsvar -psname $t_id -mode germline", true);
	}

	//annotate vcf, GSvar with frequency/depth from tumor RNA sample
	if (isset($t_rna_bam))
	{
		// annotate vcf
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency",
			"-in $variants_annotated -bam $t_rna_bam -out $variants_annotated -name rna_tum -depth", true);
		// annotate GSvar
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency",
			"-in $variants_gsvar -bam $t_rna_bam -out $variants_gsvar -name rna_tum -depth", true);
	}
}

//QCI/CGI annotation
//TODO: implementation for translocation files
$variants_qci = $full_prefix . "_var_qci.vcf.gz";				//CGI annotated vcf file
if (in_array("ci", $steps))
{
	//this step also annotates germline variants
	$variants_germline = dirname($n_bam)."/{$n_id}_var.vcf.gz";
	$variants_germline_gsvar = dirname($n_bam)."/{$n_id}.GSvar";

	// add QCI output
	$parser->execTool("Tools/converter_vcf2qci.php", "-in $variants_annotated -t_id $t_id -n_id $n_id -out $variants_qci -pass");
	
	// get cancer type from GenLab8 for diagnostic samples if not set
	if (!isset($cancer_type) && db_is_enabled("NGSD"))
	{
		$db_ngsd = DB::getInstance("NGSD");
		$t_info = get_processed_sample_info($db_ngsd, $t_id, false);
		if (!is_null($t_info) && $t_info["project_type"] == "diagnostic")
		{
			//check whether genlab credentials are available
			$db_hosts = get_path("db_host", false);
			if (!db_is_enabled("GL8"))
			{
				trigger_error("Warning: No credentials for Genlab db found. Using generic cancer type 'CANCER'!", E_USER_WARNING);
				$cancer_type = "CANCER";
			}
			else
			{
				//if no cancer_type is set try to resolve cancer_type from GENLAB
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
					trigger_error("Warning: There is no ICD10 diagnosis set in Genlab.",E_USER_WARNING);
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
				} 
				elseif($icd10_diagnosis != "" || empty($icd10_diagnosis))
				{
					trigger_error("Could not assign CGI Cancer acronym to ICD10-diagnosis '$icd10_diagnosis'. Using standard cancertype CANCER instead.",E_USER_WARNING);
					$cancer_type = "CANCER";
				}
				else
				{
					trigger_error("Could not find ICD10 diagnosis '$icd10_diagnosis' in icd10_cgi_dictionary.tsv. Aborting." ,E_USER_ERROR);
				}
			}
		}
	}
	else if (!isset($cancer_type)) //set cancer type to CANCER if not set for research samples
	{
		$cancer_type = "CANCER";
	}
	
	$parameters = "";
	if(file_exists($variants_annotated))
	{
		$parameters = $parameters . " -mutations $variants_annotated";
	}
	if(file_exists($cnvs))
	{
		$parameters = $parameters . " -cnas $cnvs";
	}
	if($cancer_type != "")
	{
		$parameters = $parameters . " -cancertype $cancer_type";
	}
	$parameters = $parameters . " -o_folder $out_folder";

	//if we know genes in target region we set parameter for this file
	$genes_in_target_region =  dirname($sys["target_file"]) . "/" .basename($sys["target_file"],".bed")."_genes.txt";
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
	$cgi_snv_result_file = $full_prefix . "_cgi_mutation_analysis.tsv";
	if(file_exists($cgi_snv_result_file) && file_exists($variants_gsvar))
	{
		$parameters = " -gsvar_in $variants_gsvar -cgi_snv_in $cgi_snv_result_file -out $variants_gsvar";
		$parser->execTool("NGS/cgi_snvs_to_gsvar.php",$parameters);
	}
	//annotate CGI cnv genes to cnv input file (which was originally created by CNVHunter)
	$cgi_cnv_result_file = $full_prefix . "_cgi_cnv_analysis.tsv";
	if(file_exists($cgi_cnv_result_file))
	{
		$parameters = " -cnv_in $cnvs -cnv_in_cgi $cgi_cnv_result_file -out $cnvs";
		$parser->execTool("NGS/cgi_annotate_cnvs.php",$parameters);
	}
	
	//parse germline variants
	$parameters = "-o_folder ".dirname($n_bam)." -is_germline -cancertype CANCER";
	if(file_exists($variants_germline) && $include_germline)
	{
		$parameters = $parameters . " -mutations $variants_germline";
		
		$result_send_data = $parser->execTool("NGS/cgi_send_data.php",$parameters,false);
		$error_code = $result_send_data[2];
		if($error_code != 0)
		{
			trigger_error("step \"ci\" for germline variants did not exit cleanly. Please check sample CGI files manually. cgi_send_data returned code ".$error_code,E_USER_WARNING);
		}
		
		$cgi_normal_snv_result_file = dirname($n_bam) . "/" . $n_id . "_cgi_mutation_analysis.tsv";
		if(file_exists($cgi_normal_snv_result_file))
		{
			$parameters = " -gsvar_in $variants_germline_gsvar -cgi_snv_in $cgi_normal_snv_result_file -out $variants_germline_gsvar";
			$parser->execTool("NGS/cgi_snvs_to_gsvar.php",$parameters,false);
		}
	}
	elseif($include_germline)
	{
		print("Could not create CGI report for germline variants because germline variant file $variants_germline does not exist.\n");
	}
}

//MSI calling
$msi_o_file = $full_prefix . "_msi.tsv";						//MSI
if (in_array("msi", $steps) && !$single_sample)
{
	$reference_genome = get_path("data_folder") . "/genomes/" .$n_sys['build'] . ".fa";
	
	//check whether file with loci exists in output folder
	//if not: intersect with loci file of reference
	$reference_loci_file = get_path("data_folder") . "/dbs/MANTIS/".$n_sys['build']."_msi_loci.bed";
	if(!file_exists($reference_loci_file))
	{
		print("Could not find loci reference file $reference_loci_file. Trying to generate it.\n");
		$parser->exec(get_path("mantis")."/tools/RepeatFinder","-i $reference_genome -o $reference_loci_file",false);
	}

	//file that contains MSI in target region -> is intersection of loci reference with target region
	$target_bed_file = $n_sys['target_file'];
	$target_loci_file = $full_prefix . "_msi_loci_target.bed";
	
	//target loci file is intersection of reference loci with target region
	if(!file_exists($target_loci_file))
	{
		$parameters = "-in ".$reference_loci_file." -in2 ".$target_bed_file ." -mode in -out ".$target_loci_file;
		$parser->exec(get_path("ngs-bits")."BedIntersect",$parameters,false);
	}

	$parameters = "-n_bam $n_bam -t_bam $t_bam -threads $threads -bed_file $target_loci_file -out " .$msi_o_file. " -build $reference_genome";
	
	if($sys['type'] == "WES") $parameters .= " -is_exome";
	
	$parser->execTool("NGS/detect_msi.php",$parameters);
}
elseif ($single_sample && in_array("msi",$steps))
{
	trigger_error("Calling microsatellite instabilities is only possible for tumor normal pairs",E_USER_NOTICE);
}

//DB import
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	
	$t_info = get_processed_sample_info($db_conn, $t_id, false);
	$n_info = get_processed_sample_info($db_conn, $n_id, false);
	
	if (is_null($t_info) || (!$single_sample && is_null($n_info)))
	{
		trigger_error("No database import since no valid processing ID (T: {$t_id}".($single_sample ? "" : " /N: {$n_id}").")", E_USER_WARNING);
	}
	else
	{
		// check sex
		$parser->execTool("NGS/db_check_gender.php", "-in $t_bam -pid $t_id");

		// import qcML files
		$log_db = dirname($t_bam)."/{$t_id}_log4_db.log";
		$qcmls = implode(" ", array_filter([
			dirname($t_bam)."/{$t_id}_stats_fastq.qcML",
			dirname($t_bam)."/{$t_id}_stats_map.qcML",
			$somaticqc
		], "file_exists"));
		$parser->execTool("NGS/db_import_qc.php", "-id $t_id -files $qcmls -force -min_depth 0 --log $log_db");

		// check tumor/normal flag
		if (!$t_info['is_tumor'])
		{
			trigger_error("Tumor sample $t_id is not flagged as tumor in NGSD!", E_USER_WARNING);
		}

		// analogous steps for normal sample, plus additional
		if (!$single_sample)
		{
			// check sex
			$parser->execTool("NGS/db_check_gender.php", "-in $n_bam -pid $n_id");

			// import qcML files
			$log_db = dirname($n_bam)."/{$n_id}_log4_db.log";
			$qcmls = implode(" ", array_filter([
				dirname($n_bam)."/{$n_id}_stats_fastq.qcML",
				dirname($n_bam)."/{$n_id}_stats_map.qcML",
				dirname($n_bam)."/{$n_id}_stats_vc.qcML"
			], "file_exists"));
			$parser->execTool("NGS/db_import_qc.php", "-id $n_id -files $qcmls -force -min_depth 0 --log $log_db");

			// check tumor/normal flag
			if ($n_info['is_tumor'])
			{
				trigger_error("Normal sample $n_id is flagged as tumor in NGSD!", E_USER_WARNING);
			}

			// update normal sample entry for tumor
			if (updateNormalSample($db_conn, $t_id, $n_id))
			{
				trigger_error("Updated normal sample ($n_id) for tumor ($t_id) in NGSD.", E_USER_NOTICE);
			}

			// import variants (not for WGS)
			if (file_exists($variants_gsvar) && $sys['type'] !== "WGS")
			{
				$parser->execTool("NGS/db_import_variants.php", "-id {$t_id}-{$n_id} -var {$variants_gsvar} -build ".$sys['build']." -force -mode somatic");
			}
		}
	}
}

?>