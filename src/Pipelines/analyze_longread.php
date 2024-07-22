<?php

/**
	@page analyze_longread
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("analyze_longread", "Complete NGS analysis pipeline for long-read data.");
$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);
//optional
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'name' is a valid processed sample name).", true);
$steps_all = array("ma", "vc", "cn", "sv", "re", "an", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, cn=copy-number analysis, sv=structural-variant analysis, re=repeat expansions calling, an=annotation, db=import into NGSD.", true, "ma,vc,sv,re,an,db");
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("skip_phasing", "Skip phasing of VCF and BAM files.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
$parser->addFlag("skip_wgs_check", "Skip the similarity check with a related short-read WGS sample.");
extract($parser->parse($argv));

// create logfile in output folder if no filepath is provided:
if (!file_exists($folder))
{
	exec2("mkdir -p $folder");
}
if ($parser->getLogFile() == "") $parser->setLogFile($folder."/analyze_".date("YmdHis").".log");

//init
$ngsbits = get_path("ngs-bits");

//determine processing system
$sys = load_system($system, $name);
$build = $sys['build'];
$genome = genome_fasta($build);

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//log server name
list($server) = exec2("hostname -f");
$user = exec('whoami');
$parser->log("Executed on server: ".implode(" ", $server)." as ".$user);



if (db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
	list($rc_id, $rc_vars_exist, $rc_cnvs_exist, $rc_svs_exist) = report_config($db, $name);
	if (in_array("vc", $steps) && $rc_vars_exist)
	{
		trigger_error("Skipping step 'vc' - Report configuration with small variants exists in NGSD!", E_USER_NOTICE);
		if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
	}
	if (in_array("cn", $steps) && $rc_cnvs_exist)
	{
		trigger_error("Skipping step 'cn' - Report configuration with CNVs exists in NGSD!", E_USER_NOTICE);
		if (($key = array_search("cn", $steps)) !== false) unset($steps[$key]);
	}
	if (in_array("sv", $steps) && $rc_svs_exist)
	{
		trigger_error("Skipping step 'sv' - Report configuration with SVs exists in NGSD!", E_USER_NOTICE);
		if (($key = array_search("sv", $steps)) !== false) unset($steps[$key]);
	}
}

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$build);
}

//output file names:
//mapping
$bam_file = $folder."/".$name.".bam";
$unmapped_bam_file = $folder."/".$name.".mod.unmapped.bam";
$lowcov_file = $folder."/".$name."_".$sys["name_short"]."_lowcov.bed";
//methylation 
$modkit_track = $folder."/".$name."_modkit_track.bed.gz";
$modkit_summary = $folder."/".$name."_modkit_summary.txt";
//variant calling
$vcf_file = $folder."/".$name."_var.vcf.gz";
$vcf_file_annotated = $folder."/".$name."_var_annotated.vcf.gz";
$var_file = $folder."/".$name.".GSvar";
$rohfile = $folder."/".$name."_rohs.tsv";
$baffile = $folder."/".$name."_bafs.igv";
$ancestry_file = $folder."/".$name."_ancestry.tsv";
$prs_file = $folder."/".$name."_prs.tsv";
$vcf_modcall = $folder."/".$name."_var_modcall.vcf";
//copy-number calling
$cnv_bin_size = get_path("cnv_bin_size_longread_wgs");
$cnv_file = $folder."/".$name."_cnvs_clincnv.tsv";
$cnv_file2 = $folder."/".$name."_cnvs_clincnv.seg";
//structural variant calling
$sv_vcf_file = $folder ."/". $name . "_var_structural_variants.vcf.gz";
$bedpe_file = substr($sv_vcf_file,0,-6)."bedpe";
//repeat expansions
$straglr_file = $folder."/".$name."_repeats.vcf";
//db import
$qc_fastq  = $folder."/".$name."_stats_fastq.qcML";
$qc_map  = $folder."/".$name."_stats_map.qcML";
$qc_vc  = $folder."/".$name."_stats_vc.qcML";
$qc_other  = $folder."/".$name."_stats_other.qcML";

$name_sample_ps = explode("_", $name, 2);
$sample_name = $name_sample_ps[0];

//mapping
if (in_array("ma", $steps))
{
	//determine input FASTQ and BAM files
	$unmapped_bam_files = glob("{$folder}/{$sample_name}_??.mod.unmapped.bam");
	$bam_files = glob("{$folder}/{$sample_name}_??.bam");
	$fastq_files = glob("{$folder}/{$sample_name}_??*.fastq.gz");

	// preference:
	// 1. mod unmapped bam
	// 2. regular bam
	// 3. fastq

	if ((count($unmapped_bam_files) === 0) && (count($bam_files) === 0) && (count($fastq_files) === 0))
	{
		trigger_error("Found no input read files in BAM or FASTQ format!", E_USER_ERROR);
	}
	
	// run mapping
	$mapping_minimap_options = [
		"-out {$bam_file}",
		"-sample {$name}",
		"-threads {$threads}",
		"-system {$system}",
		"-qc_fastq {$qc_fastq}",
		"-qc_map {$qc_map}"
	];
	if (count($unmapped_bam_files) > 0)
	{
		$mapping_minimap_options[] = "-in_bam " . implode(" ", $unmapped_bam_files);	
	}
	elseif (count($bam_files) > 0)
	{
		$mapping_minimap_options[] = "-in_bam " . implode(" ", $bam_files);
	}
	else
	{
		$mapping_minimap_options[] = "-in_fastq " . implode(" ", $fastq_files);
	}
	
	$parser->execTool("NGS/mapping_minimap.php", implode(" ", $mapping_minimap_options));

	// create methylation track
	if (contains_methylation($bam_file))
	{
		$args = array();
		$args[] = "-bam ".$bam_file;
		$args[] = "-bed ".$modkit_track;
		$args[] = "-summary ".$modkit_summary;
		$args[] = "-threads ".$threads;
		$args[] = "-build ".$build;
		$parser->execTool("NGS/vc_modkit.php", implode(" ", $args));
	}

	//low-coverage report
	$parser->exec("{$ngsbits}BedLowCoverage", "-in ".$sys['target_file']." -bam $bam_file -out $lowcov_file -cutoff 20 -threads {$threads} -ref {$genome}", true);
	if (db_is_enabled("NGSD"))
	{
		$parser->exec("{$ngsbits}BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", true);
	}

}
else
{
	//check genome build of BAM
	check_genome_build($bam_file, $build);
	
	$local_bamfile = $bam_file;
}

//variant calling
if (in_array("vc", $steps))
{
	//determine basecall model
	$basecall_model = get_basecall_model($bam_file);
	$basecall_model_path = "";

	if ($basecall_model == "")
	{
		//if no entry in BAM header -> use default model
		if (($sys["name_short"] == "SQK-LSK114") || ($sys["name_short"] == "LR-ONT-SQK-LSK114"))
		{
			$basecall_model_path = get_path("clair3_models")."/r1041_e82_400bps_hac_g632/";
		}
		else if ($sys["name_short"] == "SQK-LSK109")
		{
			$basecall_model_path = get_path("clair3_models")."/r941_prom_hac_g360+g422/";
		}
		else
		{
			trigger_error("Unsupported processing system '".$sys["shortname"]."' provided!", E_USER_ERROR);
		}

		trigger_error("No basecall info found in BAM file. Using default model at '{$basecall_model_path}'.", E_USER_NOTICE);
	}
	//special handling of old models
	else if ($basecall_model == "dna_r10.4.1_e8.2_400bps_hac@v3.5.2")
	{
		$basecall_model_path = get_path("clair3_models")."/r1041_e82_400bps_hac_g632/";

		trigger_error("Basecall info found in BAM file ('{$basecall_model}'). Using model at '{$basecall_model_path}' (special case).", E_USER_NOTICE);
	} 
	else if ($basecall_model == "dna_r9.4.1_e8_hac@v3.3")
	{
		$basecall_model_path = get_path("clair3_models")."/r941_prom_hac_g360+g422/";

		trigger_error("Basecall info found in BAM file ('{$basecall_model}'). Using model at '{$basecall_model_path}' (special case).", E_USER_NOTICE);
	} 
	//all new models can be directly derived
	else 
	{
		$reformated_model_str = strtr($basecall_model, array("dna_" => "", "." => "", "@" => "_"));
		$basecall_model_path = get_path("clair3_models")."/".$reformated_model_str."/";

		trigger_error("Basecall info found in BAM file ('{$basecall_model}'). Using model at '{$basecall_model_path}' (automatically derived).", E_USER_NOTICE);
	}

	//check if selected model is available
	if (!file_exists($basecall_model_path)) trigger_error("Basecall model at '{$basecall_model_path}' not found!", E_USER_ERROR);

	//prepare clair command
	$args = [];
	$args[] = "-bam ".$bam_file;
	$args[] = "-folder ".$folder;
	$args[] = "-name ".$name;
	$args[] = "-target ".$sys['target_file'];
	$args[] = "-target_extend 200";
	$args[] = "-threads ".$threads;
	$args[] = "-build ".$build;
	$args[] = "--log ".$parser->getLogFile();
	$args[] = "-model ".$basecall_model_path;
	
	$parser->execTool("NGS/vc_clair.php", implode(" ", $args));	

	//create b-allele frequency file
	$params = array();
	$params[] = "-vcf {$vcf_file}";
	$params[] = "-name {$name}";
	$params[] = "-out {$baffile}";
	$params[] = "-downsample 100";
	$parser->execTool("NGS/baf_germline.php", implode(" ", $params));
}

//copy-number analysis
if (in_array("cn", $steps))
{
	//create reference folder if it does not exist
	$ref_folder = get_path("data_folder")."/coverage/".$sys['name_short']."/";
	if (!is_dir($ref_folder))
	{
		mkdir($ref_folder);
		if (!chmod($ref_folder, 0777))
		{
			trigger_error("Could not change privileges of folder '{$ref_folder}'!", E_USER_ERROR);
		}
	}

	$cov_folder = $ref_folder;

	//Calling ClinCNV
	//create folder for binned coverage data - if missing
	$bin_folder = "{$ref_folder}/bins{$cnv_bin_size}/";
	if (!is_dir($bin_folder))
	{
		mkdir($bin_folder);
		if (!chmod($bin_folder, 0777))
		{
			trigger_error("Could not change privileges of folder '{$bin_folder}'!", E_USER_ERROR);
		}
	}
	$cov_folder = $bin_folder;

	
	//create BED file with GC and gene annotations - if missing
	$bed = $ref_folder."/bins{$cnv_bin_size}.bed";
	if (!file_exists($bed))
	{
		$pipeline = [
				["{$ngsbits}BedChunk", "-in ".$sys['target_file']." -n {$cnv_bin_size}"],
				["{$ngsbits}BedAnnotateGC", "-clear -ref ".$genome],
				["{$ngsbits}BedAnnotateGenes", "-out {$bed}"]
			];
		$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
	}

	//create coverage profile
	$tmp_folder = $parser->tempFolder();
	$cov_file = $cov_folder."/{$name}.cov";
	if (!file_exists($cov_file))
	{

		$parser->log("Calculating coverage file for CN calling...");
		$cov_tmp = $tmp_folder."/{$name}.cov";
		$parser->exec("{$ngsbits}BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam {$bam_file} -in {$bed} -out {$cov_tmp} -threads {$threads}", true);
		
		//copy coverage file to reference folder if valid
		if (db_is_enabled("NGSD") && is_valid_ref_sample_for_cnv_analysis($name))
		{
			$parser->log("Moving coverage file to reference folder...");
			$parser->moveFile($cov_tmp, $cov_file);
		}
		else
		{
			$cov_file = $cov_tmp;
		}
	}
	else
	{
		$parser->log("Using previously calculated coverage file for CN calling: $cov_file");
	}
	
	//perform CNV analysis
	$cnv_out = $tmp_folder."/output.tsv";
	$cnv_out2 = $tmp_folder."/output.seg";
	$args = array(
		"-cov {$cov_file}",
		"-cov_folder {$cov_folder}",
		"-bed {$bed}",
		"-out {$cnv_out}",
		"-threads {$threads}",
		"--log ".$parser->getLogFile(),
		"-max_cnvs 2000",
		"-mosaic"
	);
	
	$parser->execTool("NGS/vc_clincnv_germline.php", implode(" ", $args), true);
	
	//copy results to output folder
	if (file_exists($cnv_out)) $parser->moveFile($cnv_out, $cnv_file);
	if (file_exists($cnv_out2)) $parser->moveFile($cnv_out2, $cnv_file2);
	$mosaic = $folder."/".$name."_mosaic_cnvs.tsv";
	$sample_cnv_name = substr($cnv_out,0,-4);
	$mosaic_out = $sample_cnv_name."_mosaic.tsv";
	if (file_exists($mosaic_out)) $parser->moveFile($mosaic_out, $mosaic);

}

//structural variants
if (in_array("sv", $steps))
{
	//run Sniffles
	$parser->execTool("NGS/vc_sniffles.php", "-bam {$bam_file} -sample_ids {$name} -out {$sv_vcf_file} -threads {$threads} -build {$build}");
				
}

//phasing
if (!$skip_phasing && (in_array("vc", $steps) || in_array("sv", $steps)))
{
	//check for methylation
	$contains_methylation = contains_methylation($bam_file);

	//create modcall file
	if($contains_methylation)
	{
		//delete prev modcall file
		if(file_exists($vcf_modcall)) $parser->exec("rm", $vcf_modcall);
		$args = array();
		$args[] = "modcall";
		$args[] = "-b {$bam_file}";
		$args[] = "-r {$genome}";
		$args[] = "-t {$threads}";
		$args[] = "-o ".substr($vcf_modcall, 0, -4);
		$parser->exec(get_path("longphase"), implode(" ", $args));

		trigger_error("Methylation annotation detected. Using intermediate modcall step.", E_USER_NOTICE);
	}

	//run phasing by LongPhase on VCF files
	$phased_tmp = $parser->tempFile(".vcf", "longphase");
	$phased_sv_tmp = substr($phased_tmp,0,-4)."_SV.vcf";
	$args = array();
	$args[] = "phase";
	$args[] = "-s {$vcf_file}";
	$args[] = "-b {$bam_file}";
	$args[] = "-r {$genome}";
	$args[] = "-t {$threads}";
	$args[] = "-o ".substr($phased_tmp, 0, -4);
	$args[] = "--ont";
	$args[] = "--indels";
	if (file_exists($sv_vcf_file)) $args[] = "--sv-file {$sv_vcf_file}";
	if ($contains_methylation) $args[] = "--mod-file {$vcf_modcall}";
	

	$parser->exec(get_path("longphase"), implode(" ", $args));
	
	//create compressed file and index
	$parser->exec("bgzip", "-c $phased_tmp > {$vcf_file}", false);
	$parser->exec("tabix", "-f -p vcf $vcf_file", false);
	if (file_exists($sv_vcf_file))
	{
		$parser->exec("bgzip", "-c $phased_sv_tmp > {$sv_vcf_file}", false);
		$parser->exec("tabix", "-f -p vcf $sv_vcf_file", false);
	}

	//tag BAM file 
	$args = array();
	$tagged_bam_file = $parser->tempFile(".tagged.bam");
	$args[] = "haplotag";
	$args[] = "-s {$vcf_file}";
	$args[] = "-b {$bam_file}";
	$args[] = "-r {$genome}";
	$args[] = "-t {$threads}";
	if (file_exists($sv_vcf_file)) $args[] = "--sv-file {$sv_vcf_file}";
	$args[] = "-o ".substr($tagged_bam_file, 0, -4);

	$parser->exec(get_path("longphase"), implode(" ", $args));
	$parser->indexBam($tagged_bam_file, $threads);

	//TODO: compare tagged BAM with input BAM

	//replace current bam file
	$parser->copyFile($tagged_bam_file, $bam_file);
	$parser->copyFile($tagged_bam_file.".bai", $bam_file.".bai");

	// check if copy was successful
	if (!file_exists($bam_file) || filesize($tagged_bam_file) != filesize($bam_file))
	{
		trigger_error("Error during coping BAM file! File sizes don't match!", E_USER_ERROR);
	}
}

// repeat expansion
if (in_array("re", $steps))
{
	//Repeat-expansion calling using straglr
	$variant_catalog = repository_basedir()."/data/repeat_expansions/straglr_variant_catalog_grch38.bed";
	$parser->execTool("NGS/vc_straglr.php", "-in {$bam_file} -out {$straglr_file} -loci {$variant_catalog} -threads {$threads} -build {$build}");
}

// annotation
if (in_array("an", $steps))
{
	//annotate small variant VCF file
	if (file_exists($vcf_file))
	{
		check_genome_build($vcf_file, $build);

		//annotation
		$args = [];
		$args[] = "-out_name ".$name;
		$args[] = "-out_folder ".$folder;
		$args[] = "-system ".$system;
		$args[] = "--log ".$parser->getLogFile();
		$args[] = "-threads ".$threads;

		//run annotation pipeline
		$parser->execTool("Pipelines/annotate.php", implode(" ", $args));

		//ROH detection
		$args = [];
		$args[] = "-in $vcf_file_annotated";
		$args[] = "-out $rohfile";
		$args[] = "-var_af_keys gnomADg_AF";
		$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; //optional because of license
		$args[] = "-annotate ".repository_basedir()."/data/gene_lists/genes.bed ".(file_exists($omim_file) ? $omim_file : "");
		$parser->exec("{$ngsbits}RohHunter", implode(" ", $args), true);

		//PRS calculation 
		$prs_folder = repository_basedir()."/data/misc/prs/";
		$prs_scoring_files = glob($prs_folder."/*_".$build.".vcf");
		if (count($prs_scoring_files) > 0)
		{
			$parser->exec("{$ngsbits}VcfCalculatePRS", "-in {$vcf_file} -bam {$bam_file} -out $prs_file -prs ".implode(" ", $prs_scoring_files)." -ref $genome", true);
		}

		//determine ancestry
		if (ngsbits_build($build) != "non_human")
		{
			$parser->exec($ngsbits."SampleAncestry", "-in {$vcf_file} -out {$ancestry_file} -build ".ngsbits_build($build), true);
		}
	}

	// annotate CNV file
	if (file_exists($cnv_file))
	{
		$repository_basedir = repository_basedir();
		$data_folder = get_path("data_folder");
		$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnv_file} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$cnv_file}", true);
		$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnv_file} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$cnv_file}", true);
		$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnv_file} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed -no_duplicates -url_decode -out {$cnv_file}", true);
		$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnv_file} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2023-07.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$cnv_file}", true);

		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2023_3.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnv_file} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$cnv_file}", true);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnv_file} -in2 {$omim_file} -no_duplicates -url_decode -out {$cnv_file}", true);
		}

		//annotate additional gene info
		$parser->exec($ngsbits."CnvGeneAnnotation", "-in {$cnv_file} -add_simple_gene_names -out {$cnv_file}", true);
		// skip annotation if no connection to the NGSD is possible
		if (db_is_enabled("NGSD"))
		{
			//annotate overlap with pathogenic CNVs
			$parser->exec($ngsbits."NGSDAnnotateCNV", "-in {$cnv_file} -out {$cnv_file}", true);
		}
	}
	else
	{
		trigger_error("CNV file {$cnv_file} does not exist, skipping CNV annotation!", E_USER_WARNING);
	}

	//annotate SV VCF file
	if (file_exists($sv_vcf_file))
	{
		check_genome_build($sv_vcf_file, $build);

		//create BEDPE files
		$parser->exec("{$ngsbits}VcfToBedpe", "-in $sv_vcf_file -out $bedpe_file", true);

		//add gene info annotation
		if (db_is_enabled("NGSD"))
		{
			$parser->exec("{$ngsbits}BedpeGeneAnnotation", "-in $bedpe_file -out $bedpe_file -add_simple_gene_names", true);
		}

		//add NGSD counts from flat file
		$ngsd_annotation_folder = get_path("data_folder")."/dbs/NGSD/";
		$ngsd_sv_files = array("sv_deletion.bedpe.gz", "sv_duplication.bedpe.gz", "sv_insertion.bedpe.gz", "sv_inversion.bedpe.gz", "sv_translocation.bedpe.gz");
		$db_file_dates = array();

		// check file existance
		$all_files_available = file_exists($ngsd_annotation_folder."sv_breakpoint_density.igv");
		foreach ($ngsd_sv_files as $filename) 
		{
			if(!(file_exists($ngsd_annotation_folder.$filename) && file_exists($ngsd_annotation_folder.$filename.".tbi")))
			{
				$all_files_available = false;
				break;
			}
		}
		if ($all_files_available)
		{
			// store flat file modification date to detect changes during annotation 
			foreach ($ngsd_sv_files as $filename)
			{
				$db_file_dates[$filename] = filemtime($ngsd_annotation_folder.$filename);
				if ($db_file_dates[$filename] == false)
				{
					trigger_error("Cannot get modification date of '".$ngsd_annotation_folder.$filename."'!",E_USER_ERROR);
				}
			}
			
			//perform annotation
			$parser->exec("{$ngsbits}BedpeAnnotateCounts", "-in $bedpe_file -out $bedpe_file -processing_system ".$sys["name_short"]." -ann_folder {$ngsd_annotation_folder}", true);
			$sys_specific_density_file = $ngsd_annotation_folder."sv_breakpoint_density_".$sys["name_short"].".igv";
			if (file_exists($sys_specific_density_file))
			{
				$parser->exec("{$ngsbits}BedpeAnnotateBreakpointDensity", "-in {$bedpe_file} -out {$bedpe_file} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv -density_sys {$sys_specific_density_file}", true);
			}
			else
			{
				$parser->exec("{$ngsbits}BedpeAnnotateBreakpointDensity", "-in {$bedpe_file} -out {$bedpe_file} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv", true);
			}
			

			// check if files changed during annotation
			foreach ($ngsd_sv_files as $filename)
			{
				if ($db_file_dates[$filename] != filemtime($ngsd_annotation_folder.$filename))
				{
					trigger_error("Annotation file '".$ngsd_annotation_folder.$filename."' has changed during annotation!",E_USER_ERROR);
				}
			}

		}
		else
		{
			trigger_error("Cannot annotate NGSD counts! At least one required file in '{$ngsd_annotation_folder}' is missing!", E_USER_WARNING);
		}
		
		//add optional OMIM annotation
		$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; 
		if(file_exists($omim_file))//OMIM annotation (optional because of license)
		{
			$parser->exec("{$ngsbits}BedpeAnnotateFromBed", "-in $bedpe_file -out $bedpe_file -bed $omim_file -url_decode -replace_underscore -col_name OMIM", true);
		}

		//add CNV overlap annotation
		if (file_exists($cnv_file))
		{
			$parser->exec("{$ngsbits}BedpeAnnotateCnvOverlap", "-in $bedpe_file -out $bedpe_file -cnv $cnv_file", true);
		}

		//write genotype in own column
		$parser->exec("{$ngsbits}BedpeExtractGenotype", "-in $bedpe_file -out $bedpe_file -include_unphased", true);

		//extract columns
		$parser->exec("{$ngsbits}BedpeExtractInfoField", "-in $bedpe_file -out $bedpe_file -info_fields SVLEN,SUPPORT,COVERAGE,AF", true);

		//update sample entry 
		update_gsvar_sample_header($bedpe_file, array($name=>"Affected"));		
	}

}

// collect other QC terms - if CNV or SV calling was done
if ((in_array("cn", $steps) || in_array("sv", $steps) || in_array("an", $steps)))
{
	$terms = [];
	$sources = [];
	
	//CNVs
	if (file_exists($cnv_file))
	{
		$cnv_count_hq = 0;
		$cnv_count_hq_autosomes = 0;
		$cnv_count_loss = 0;
		$cnv_count_gain = 0;
		$h = fopen2($cnv_file, 'r');
		while(!feof($h))
		{
			$line = trim(fgets($h));
			if ($line=="") continue;
			
			if (starts_with($line, "##mean correlation to reference samples:"))
			{
				$value = trim(explode(":", $line)[1]);
				$terms[] = "QC:2000114\t{$value}";
			}
			if (starts_with($line, "##number of iterations:"))
			{
				$value = trim(explode(":", $line)[1]);
				$terms[] = "QC:2000115\t{$value}";
			}
			
			if ($line[0]!="#")
			{
				$parts = explode("\t", $line);
				$ll = $parts[4];
				if ($ll>=20)
				{
					++$cnv_count_hq;
					
					$chr = $parts[0];
					if (is_numeric(strtr($chr, ["chr"=>""])))
					{
						++$cnv_count_hq_autosomes;
						$cn = $parts[3];
						if ($cn<2) ++$cnv_count_loss;
						if ($cn>2) ++$cnv_count_gain;
					}
				}
			}
		}
		fclose($h);
		
		//counts (all, loss, gain)
		$terms[] = "QC:2000113\t{$cnv_count_hq}";
		if ($cnv_count_hq_autosomes>0)
		{
			$terms[] = "QC:2000118\t".number_format(100.0*$cnv_count_loss/$cnv_count_hq_autosomes, 2);
			$terms[] = "QC:2000119\t".number_format(100.0*$cnv_count_gain/$cnv_count_hq_autosomes, 2);
		}
		
		$sources[] = $cnv_file;
	}

	//SVs
	if (file_exists($bedpe_file))
	{
		$sv_count_pass = 0;
		$sv_count_del = 0;
		$sv_count_dup = 0;
		$sv_count_ins = 0;
		$sv_count_inv = 0;
		$sv_count_bnd = 0;
		$h = fopen2($bedpe_file, 'r');
		while(!feof($h))
		{
			$line = trim(fgets($h));
			if ($line=="" || $line[0]=="#") continue;
			
			
			$parts = explode("\t", $line);
			$filter = trim($parts[11]);
			if ($filter=="PASS")
			{
				++$sv_count_pass;
				$type = trim($parts[10]);
				if ($type=="DEL") ++$sv_count_del;
				if ($type=="DUP") ++$sv_count_dup;
				if ($type=="INS") ++$sv_count_ins;
				if ($type=="INV") ++$sv_count_inv;
				if ($type=="BND") ++$sv_count_bnd;
				
			}
		}
		fclose($h);
		
		$terms[] = "QC:2000117\t{$sv_count_pass}";
		if ($sv_count_pass>0)
		{
			$terms[] = "QC:2000120\t".number_format(100.0*$sv_count_del/$sv_count_pass, 2);
			$terms[] = "QC:2000121\t".number_format(100.0*$sv_count_dup/$sv_count_pass, 2);
			$terms[] = "QC:2000122\t".number_format(100.0*$sv_count_ins/$sv_count_pass, 2);
			$terms[] = "QC:2000123\t".number_format(100.0*$sv_count_inv/$sv_count_pass, 2);
			$terms[] = "QC:2000124\t".number_format(100.0*$sv_count_bnd/$sv_count_pass, 2);
		}
		
		$sources[] = $bedpe_file;
	}
	
	if(count($sources) > 0)
	{
		//create qcML file
		$tmp = $parser->tempFile("qc.tsv");
		file_put_contents($tmp, implode("\n", $terms));
		$parser->exec("{$ngsbits}TsvToQC", "-in $tmp -out $qc_other -sources ".implode(" ", $sources));
	}
}

// create Circos plot - if small variant, CNV or SV calling was done
if (in_array("an", $steps))
{
	if (file_exists($cnv_file))
	{
		if (file_exists($cnv_file2))
		{
			$parser->execTool("NGS/create_circos_plot.php", "-folder $folder -name $name -build ".$build);
		}
		else
		{
			trigger_error("CNV file $cnv_file2 missing. Cannot create Circos plot!", E_USER_WARNING);
		}
	}
	else
	{
		trigger_error("CNV file $cnv_file missing. Cannot create Circos plot!", E_USER_WARNING);
	}
}

//import to database
if (in_array("db", $steps))
{
	//import ancestry
	if(file_exists($ancestry_file)) $parser->execTool("NGS/db_import_ancestry.php", "-id {$name} -in {$ancestry_file}");
	
	//import QC
	$qc_files = array($qc_fastq, $qc_map);
	if (file_exists($qc_vc)) $qc_files[] = $qc_vc; 
	if (file_exists($qc_other)) $qc_files[] = $qc_other; 
	$parser->exec("{$ngsbits}/NGSDImportSampleQC", "-ps $name -files ".implode(" ", $qc_files)." -force");
	
	//check gender
	$parser->execTool("NGS/db_check_gender.php", "-in $bam_file -pid $name");	
	

	if (!$skip_wgs_check)
	{
		//check similarity to (sr) WGS sample
		list($stdout, $stderr, $exit_code) = $parser->exec("{$ngsbits}/NGSDSameSample", "-ps $name -system_type WGS");
		//parse stdout
		$same_samples = array();
		foreach ($stdout as $line) 
		{
			if (starts_with($line, "#")) continue;
			if (trim($line) == "") continue;
			$same_samples[] = trim(explode("\t", $line)[0]);
		}
		if (count($same_samples) < 1)
		{
			trigger_error("No related WGS sample found for {$name}. Cannot perform similarity check!", E_USER_WARNING);
		} 
		else
		{
			//check same-sample correlation for each WGS sample
			$min_corr = 0.85;

			foreach ($same_samples as $processed_sample) 
			{
				$ps_gsvar = trim($parser->exec(get_path("ngs-bits")."SamplePath", "-ps {$processed_sample} -type GSVAR", true)[0][0]);
				if (!file_exists($ps_gsvar))
				{
					trigger_error("GSvar file {$ps_gsvar} not found! Skipping sample similarity check", E_USER_WARNING);
					continue;
				}
				$output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", "-in {$var_file} {$ps_gsvar} -mode gsvar -build ".ngsbits_build($sys['build']), true);
				$correlation = explode("\t", $output[0][1])[3];
				if ($correlation<$min_corr)
				{
					trigger_error("The genotype correlation of {$name} and {$processed_sample} is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
				}
				else
				{
					trigger_error("The genotype correlation of {$name} and {$processed_sample} is {$correlation}.", E_USER_NOTICE);
				}
			}
		}
	}
	

	//import variants
	$args = ["-ps {$name}"];
	$import = false;
	if (file_exists($var_file))
	{
		//check genome build
		check_genome_build($var_file, $build);
		
		$args[] = "-var {$var_file}";
		$args[] = "-var_force";
		$import = true;
	}
	if (file_exists($cnv_file))
	{
		//check genome build
		//this is not possible for CNVs because the file does not contain any information about it
		
		$args[] = "-cnv {$cnv_file}";
		$args[] = "-cnv_force";
		$import = true;
	}
	if (file_exists($bedpe_file))
	{
		//check genome build
		check_genome_build($bedpe_file, $build);
		
		$args[] = "-sv {$bedpe_file}";
		$args[] = "-sv_force";
		$import = true;
	}
	if (file_exists($straglr_file))
	{
		//check genome build
		check_genome_build($straglr_file, $build);
		
		$args[] = "-re {$straglr_file}";
		$args[] = "-re_force";
		$import = true;
	}
	if ($import)
	{
		$parser->exec("{$ngsbits}NGSDAddVariantsGermline", implode(" ", $args), true);
	}
}

?>
