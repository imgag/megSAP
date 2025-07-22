<?php

/**
 * @page somatic_tumor_only
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_tumor_only", "Tumor-only analysis pipeline.");

//mandatory
$parser->addInfile("t_bam", "Tumor sample BAM file.", false);
$parser->addString("out_folder", "Output folder.", false);

//optional
$parser->addInfile("t_rna_bam", "Tumor RNA sample BAM file.", true);
$parser->addString("prefix", "Output file prefix.", true, "somatic");

$steps_all = array("vc", "vi", "cn", "an", "an_rna", "db", "msi");
$parser->addString("steps", "Comma-separated list of steps to perform:\n" .
	"vc=variant calling, an=annotation,\n" .
	"cn=copy-number analysis\n".
	"an_rna=annotate data from somatic RNA files,\n".
	"vi=virus detection, db=database import",
	true, "vc,cn,an,vi,an_rna,db");

$parser->addInfile("system", "Processing system file used for tumor DNA sample (resolved from NGSD via tumor BAM by default).", true);
$parser->addFlag("skip_contamination_check", "Skips check of female tumor sample for male SRY DNA.");
$parser->addFlag("skip_correlation", "Skip sample correlation check.");
$parser->addFlag("skip_low_cov", "Skip low coverage statistics.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");

//default cut-offs
$parser->addFloat("min_af", "Allele frequency detection limit (for tumor-only calling only).", true, 0.05);
$parser->addFloat("min_depth_t", "Tumor sample coverage cut-off for low coverage statistics.", true, 60);
$parser->addString("rna_ref_tissue", "Reference data for RNA annotation", true);
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

if (in_array("msi",$steps))
{
	trigger_error("Calling microsatellite instabilities is only possible for tumor normal pairs",E_USER_NOTICE);
}


###################################### SCRIPT START ######################################
if (!file_exists($out_folder))
{
	exec2("mkdir -p $out_folder");
}

$out_folder = realpath($out_folder);

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($out_folder."/somatic_tumor_only".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//output prefix
$full_prefix = "{$out_folder}/{$prefix}";

//IDs, system and target region
$t_id = basename2($t_bam);
$t_basename = dirname($t_bam)."/".$t_id;
$sys = load_system($system, $t_id);
$roi = $sys["target_file"];
$ref_genome = genome_fasta($sys['build']);

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);
}

check_genome_build($t_bam, $sys['build']);

if (!$skip_contamination_check)
{
	$parser->execTool("Auxilary/check_sry_coverage.php", "-t_bam $t_bam -build ".$sys['build']);
}
else trigger_error("Skipping check of female tumor sample $t_bam for contamination with male genomic DNA.", E_USER_WARNING);

//low coverage statistics
$low_cov = "{$full_prefix}_stat_lowcov.bed";					// low coverage BED file
if ($sys['type'] !== "WGS" && !empty($roi) && !$skip_low_cov)
{
	$parser->execApptainer("ngs-bits", "BedLowCoverage", "-in ".realpath($roi)." -bam $t_bam -out $low_cov -cutoff $min_depth_t -threads {$threads} -ref {$ref_genome}", [$roi, $t_bam, $ref_genome], [dirname($low_cov)]);

	// annotate with gene names
	if (db_is_enabled("NGSD"))
	{
		$parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-in $low_cov -extend 25 -out $low_cov", [$low_cov]);
	}
}

//variant calling
$manta_indels  = $full_prefix . "_manta_var_smallIndels.vcf.gz";	// small indels from manta
$manta_sv      = $full_prefix . "_manta_var_structural.vcf.gz";		// structural variants (vcf)
$manta_sv_bedpe= $full_prefix . "_manta_var_structural.bedpe"; 		// structural variants (bedpe)
$variants      = $full_prefix . "_var.vcf.gz";					// variants
$ballele       = $full_prefix . "_bafs.igv";					// B-allele frequencies
$hla_file_tumor = "{$t_basename}_hla_genotyper.tsv";

if (in_array("vc", $steps))
{
	//get HLA types
	$parser->execTool("Tools/hla_genotyper.php", "-bam $t_bam -name $t_id -out ".$hla_file_tumor);
	
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
			"-t_bam {$t_bam}",
			"-out {$manta_sv}",
			"-build ".$sys['build'],
			"-smallIndels {$manta_indels}",
			"-threads {$threads}"
		];
		if ($sys['type'] !== "WGS") //use exome flag for non targeted / exome samples (i.e. non WGS samples)
		{
			$args_manta[] = "-exome";
		}
		if (!empty($roi))
		{
			$args_manta[] = "-target {$roi}";
		}
		$parser->execTool("Tools/vc_manta.php", implode(" ", $args_manta));
		
		$parser->execApptainer("ngs-bits", "VcfToBedpe", "-in $manta_sv -out $manta_sv_bedpe", [$manta_sv], [dirname($manta_sv_bedpe)]);
		
		if( db_is_enabled("NGSD") )
		{
			$parser->execApptainer("ngs-bits", "BedpeGeneAnnotation", "-in $manta_sv_bedpe -out $manta_sv_bedpe -add_simple_gene_names", [$manta_sv_bedpe]);
		}
	}
	
	//Varscan2 calling
	$parser->execTool("Tools/vc_varscan2.php", "-bam $t_bam -out $variants -build " .$sys['build']. " -target ". $roi. " -name $t_id -min_af $min_af");
	$parser->execApptainer("htslib", "tabix", "-f -p vcf $variants", [], [dirname($variants)]);

	//add somatic BAF file
	//create b-allele frequency file
	$params = array();
	$params[] = "-vcf {$variants}";
	$params[] = "-name {$t_id}";
	$params[] = "-out {$ballele}";
	
	if ($sys['type'] === "WGS")
	{
		$baf_args[] = "-downsample 100";
	}
	
	$parser->execTool("Tools/baf_germline.php", implode(" ", $params));
}

//Viral sequences alignment
$viral         = "{$t_basename}_viral.tsv";					// viral sequences results
$viral_bam     = "{$t_basename}_viral.bam";					// viral sequences alignment
$viral_bam_raw = "{$t_basename}_viral_before_dedup.bam";		// viral sequences alignment (no deduplication)
$viral_bed     = get_path("data_folder") . "/enrichment/somatic_viral.bed"; //viral enrichment
$viral_genome  = get_path("data_folder") . "/genomes/somatic_viral.fa"; //viral reference genome
if (in_array("vi", $steps))
{
	if(!file_exists($viral_genome) || !file_exists($viral_bed))
	{
		trigger_error("Could not find viral reference genome {$viral_genome} or target file {$viral_bed}. Skipping step \"vi\".", E_USER_WARNING);
	}
	else
	{
		//detection of viral sequences
		$t_bam_dedup = "{$t_basename}_before_dedup.bam";
		$t_bam_map_qc = "{$t_basename}_stats_map.qcML";
		$dedup_used = file_exists($t_bam_dedup);
		$vc_viral_args = [
			"-in ".($dedup_used ? $t_bam_dedup : $t_bam),
			"-viral_bam ".$viral_bam,
			"-viral_bam_raw ".$viral_bam_raw,
			"-viral_cov ".$viral,
			"-viral_chrs chrNC_007605",
			"-build_viral somatic_viral",
			"-in_qcml ".$t_bam_map_qc,
			"-threads ".$threads
		];
		if ($dedup_used) $vc_viral_args[] = "-barcode_correction";
		if ($no_sync) $vc_viral_args[] = "no_sync";
		$parser->execTool("Tools/vc_viral_load.php", implode(" ", $vc_viral_args));
	}
}

//CNV calling
$som_clincnv = $full_prefix . "_clincnv.tsv"; //ClinCNV output file
if(in_array("cn",$steps))
{
	// copy number variant calling
	$tmp_folder = $parser->tempFolder();
	//directory with reference coverage files and directory with reference off-target coverage files
	$ref_folder_t = get_path("data_folder")."/coverage/".$sys['name_short']."-tumor";
	$ref_folder_n = get_path("data_folder")."/coverage/".$sys['name_short'];
	$ref_folder_t_off_target = $ref_folder_t . "_off_target";
	
	create_directory($ref_folder_t);
	create_directory($ref_folder_t_off_target);
	
	//Generate file with off-target region
	$off_target_bed = get_path("data_folder")."/coverage/off_target_beds/".$sys['name_short'].".bed";
	if(!file_exists($off_target_bed)) create_off_target_bed_file($off_target_bed,$sys['target_file'], $ref_genome);
	
	$target_bed = "";
	
	//get correct coverage files, target bed and off-target bed based on system type:
	if ($sys["type"] != "WGS" && $sys['type'] != "WGS (shallow)")
	{
		$target_bed = $ref_folder_t."/roi_annotated.bed";
		if (!file_exists($target_bed))
		{
			$pipeline = [
					["", $parser->execApptainer("ngs-bits", "BedAnnotateGC", "-in ".realpath($sys['target_file'])." -clear -ref {$ref_genome}", [$sys['target_file'], $ref_genome], [], true)],
					["", $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-out {$target_bed}", [], [dirname($target_bed)], true)],
				];
			$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
		}
	}
	else // WGS sample
	{
		$bin_size = get_path("cnv_bin_size_wgs");
		
		//target bed
		$target_bed = $ref_folder_t."/bins{$bin_size}.bed"; // is outside of bins folder
		if (!file_exists($target_bed))
		{
			$pipeline = [
					["", $parser->execApptainer("ngs-bits", "BedChunk", "-in ".realpath($sys['target_file'])." -n {$bin_size}", [$sys['target_file']], [], true)],
					["", $parser->execApptainer("ngs-bits", "BedAnnotateGC", "-clear -ref {$ref_genome}", [$ref_genome], [], true)],
					["", $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-out {$target_bed}", [], [dirname($target_bed)], true)]
				];
			$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
		}
		
		//coverage files:
		
		$bin_folder_t = "{$ref_folder_t}/bins{$bin_size}/";
		$bin_folder_n = "{$ref_folder_n}/bins{$bin_size}/";
		create_directory($bin_folder_t);
		create_directory($bin_folder_n);
		$ref_folder_n = $bin_folder_n;
		$ref_folder_t = $bin_folder_t;		
	}
	
	
	/***************************************************
	 * GENERATE AND COPY COVERAGE FILES TO DATA FOLDER *
	 ***************************************************/
	// coverage for tumor sample + off-target files
	$t_cov = "{$tmp_folder}/{$t_id}.cov";
	$t_cov_off_target = "{$tmp_folder}/{$t_id}_off_target.cov";
	$ref_file_t = "{$ref_folder_t}/{$t_id}.cov";
	$ref_file_t_off_target = "{$ref_folder_t_off_target}/{$t_id}.cov";
	
	$parser->execApptainer("ngs-bits", "BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam $t_bam -in $target_bed -out $t_cov -threads {$threads} -ref {$ref_genome}", [$t_bam, $target_bed, $ref_genome]);
	$parser->execApptainer("ngs-bits", "BedCoverage", "-clear -min_mapq 10 -decimals 4 -in $off_target_bed -bam $t_bam -out $t_cov_off_target -threads {$threads} -ref {$ref_genome}", [$off_target_bed, $t_bam, $ref_genome]);

	
	//copy tumor sample coverage file to reference folder (has to be done before ClinCNV call to avoid analyzing the same sample twice)
	if (db_is_enabled("NGSD") && is_valid_ref_tumor_sample_for_cnv_analysis($t_id))
	{
		$parser->log("Moving tumor sample coverage file to reference folder...");
		$parser->copyFile($t_cov, $ref_file_t); 
		$t_cov = $ref_file_t;

		$parser->copyFile($t_cov_off_target,$ref_file_t_off_target);
		$t_cov_off_target = $ref_file_t_off_target;
	}
	
	//perform CNV analysis
	$baf_folder = get_path("data_folder")."/coverage/". $sys['name_short']."_bafs";
	
	$args = array(
		"-cov {$t_cov}",
		"-cov_folder {$ref_folder_n}",
		"-bed {$target_bed}",
		"-out {$som_clincnv}",
		"-tumor_only",
		"-max_cnvs 200",
	);
	
	if ($sys["type"] != "WGS" && $sys["type"] != "WGS (shallow)")
	{
		$args[] = "-bed_off {$off_target_bed}";
		$args[] = "-cov_off {$t_cov_off_target}";
		$args[] = "-cov_folder_off {$ref_folder_t_off_target}";
	}
	
	if(is_dir($baf_folder))
	{
		$args[] = "-baf_folder {$baf_folder}";
	}
	$args[] = "--log ".$parser->getLogFile();
	$parser->execTool("Tools/vc_clincnv_germline.php", implode(" ", $args), true);

	// annotate CNV file
	$repository_basedir = repository_basedir();
	$data_folder = get_path("data_folder");
	$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$som_clincnv}", [$som_clincnv, "{$repository_basedir}/data/misc/cn_pathogenic.bed"]);
	$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed -no_duplicates -url_decode -out {$som_clincnv}", [$som_clincnv, "{$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed"]);
	$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2025-02.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$som_clincnv}", [$som_clincnv, "{$data_folder}/dbs/ClinVar/clinvar_cnvs_2025-02.bed"]);

	$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2024_4.bed"; //optional because of license
	if (file_exists($hgmd_file))
	{
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$som_clincnv}", [$som_clincnv, $hgmd_file]);
	}
	$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
	if (file_exists($omim_file))
	{
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$omim_file} -no_duplicates -url_decode -out {$som_clincnv}", [$som_clincnv, $omim_file]);
	}
	$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$repository_basedir}/data/gene_lists/genes.bed -no_duplicates -url_decode -out {$som_clincnv}", [$som_clincnv, "{$repository_basedir}/data/gene_lists/genes.bed"]);

	//annotate additional gene info
	if(db_is_enabled("NGSD"))
	{
		$parser->execApptainer("ngs-bits", "CnvGeneAnnotation", "-in {$som_clincnv}  -add_simple_gene_names -out {$som_clincnv}", [$som_clincnv]);
	}
	
	//annotate overlap with pathogenic CNVs
	if (db_is_enabled("NGSD"))
	{
		$parser->execApptainer("ngs-bits", "NGSDAnnotateCNV", "-in {$som_clincnv} -out {$som_clincnv}", [$som_clincnv]);
	}
}

//annotation
$variants_annotated = $full_prefix . "_var_annotated.vcf.gz";	// annotated variants
$variants_gsvar     = $full_prefix . ".GSvar";					// GSvar variants
$somaticqc          = $full_prefix . "_stats_som.qcML";			// SomaticQC qcML
$cfdna_folder       = $full_prefix . "_cfDNA_candidates";       // folder containing cfDNA monitoring variants 
if (in_array("an", $steps))
{
	check_genome_build($variants, $sys['build']);
	
	// annotate vcf (in temp folder)
	$tmp_folder1 = $parser->tempFolder();
	$tmp_vcf = "{$tmp_folder1}/{$prefix}_var_annotated.vcf.gz";
	$parser->execTool("Tools/annotate.php", "-out_name $prefix -out_folder $tmp_folder1 -system $system -vcf $variants -somatic -threads $threads");

	//add sample info to VCF header
	$s = Matrix::fromTSV($tmp_vcf);
	$comments = $s->getComments();
	$comments[] = gsvar_sample_header($t_id, array("IsTumor" => "yes"), "#", "");

	$s->setComments(sort_vcf_comments($comments));
	$s->toTSV($tmp_vcf);

	// zip and index vcf file
	$parser->execApptainer("htslib", "bgzip", "-c $tmp_vcf > $variants_annotated", [], [dirname($variants_annotated)]);
	$parser->execApptainer("htslib", "tabix", "-f -p vcf $variants_annotated", [], [dirname($variants_annotated)]);

	// convert vcf to GSvar
	$args = array("-in $tmp_vcf", "-out $variants_gsvar", "-t_col $t_id");
	$parser->execTool("Tools/vcf2gsvar_somatic.php", implode(" ", $args));
	
	//Annotate data from network of cancer genes
	$parser->execTool("Tools/an_somatic_gsvar.php" , "-gsvar_in $variants_gsvar -out $variants_gsvar -include_ncg");
}

if (in_array("an_rna", $steps))
{
	$args = array();
	$args[] = "-t_bam $t_bam";
	$args[] = "-full_prefix $full_prefix";
	$args[] = "-steps ".count($steps);
	if (file_exists($system)) $args[] = "-system $system";
	if (file_exists($t_rna_bam)) $args[] = "-t_rna_bam $t_rna_bam";
	if ($rna_ref_tissue != "")  $args[] = "-rna_ref_tissue $rna_ref_tissue";
	if ($skip_correlation) $args[] = "-skip_correlation";

	$parser->execTool("Tools/an_somatic_rna.php", implode(" ", $args));
}

//Collect QC terms if necessary
$qc_other = $full_prefix."_stats_other.qcML";
if (in_array("vc", $steps) || in_array("vi", $steps) || in_array("cn", $steps) || in_array("db", $steps))
{
	$parser->execTool("Auxilary/create_qcml.php", "-full_prefix $full_prefix -t_basename $t_basename -tumor_only");
}

//DB import
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	
	$t_info = get_processed_sample_info($db_conn, $t_id, false);
		
	if (is_null($t_info))
	{
		trigger_error("No database import since no valid processing ID (T: {$t_id})", E_USER_WARNING);
	}
	else
	{
		// import qcML files
		$qcmls = implode(" ", array_filter([
			"{$t_basename}_stats_fastq.qcML",
			"{$t_basename}_stats_map.qcML",
			$somaticqc,
			$qc_other
		], "file_exists"));
		$parser->execApptainer("ngs-bits", "NGSDImportSampleQC", "-ps $t_id -files $qcmls -force", ["{$t_basename}_stats_fastq.qcML", "{$t_basename}_stats_map.qcML", $somaticqc, $qc_other]);

		// check tumor/normal flag
		if (!$t_info['is_tumor'])
		{
			trigger_error("Tumor sample $t_id is not flagged as tumor in NGSD!", E_USER_WARNING);
		}

		// import variants into NGSD
		if (file_exists($variants_gsvar))
		{
			check_genome_build($variants_gsvar, $sys['build']);
			$parser->execApptainer("ngs-bits", "NGSDAddVariantsSomatic", "-t_ps $t_id -var $variants_gsvar -force", [$variants_gsvar]);
		}
				
		//add secondary analysis (if missing)
		$parser->execTool("Tools/db_import_secondary_analysis.php", "-type 'somatic' -gsvar {$variants_gsvar}");
	}
}
?>
