<?php

/**
	@page analyze_longread
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("analyze_longread", "Complete NGS analysis pipeline for long-read data.");
$parser->addInfile("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);
//optional
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'name' is a valid processed sample name).", true);
$steps_all = array("ma", "vc", "cn", "sv", "ph", "re", "pg", "an", "me", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, cn=copy-number analysis, sv=structural-variant analysis, ph=phasing, re=repeat expansions calling, pg=paralogous genes calling, an=annotation, me=methylation calling, db=import into NGSD.", true, implode(",", $steps_all));
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
$parser->addFlag("no_gender_check", "Skip gender check (done between mapping and variant calling).");
$parser->addFlag("skip_wgs_check", "Skip the similarity check with a related short-read WGS sample.");
$parser->addFlag("bam_output", "Output file format for mapping is BAM instead of CRAM.");
$parser->addFlag("gpu", "Use GPU supported tools where possible (currently DeepVariant and Clair3 in step 'vc').");
$parser->addFloat("min_af", "Minimum VAF cutoff used for variant calling (DeepVariant 'min-alternate-fraction' parameter).", true, 0.1);
$parser->addFloat("min_bq", "Minimum base quality used for variant calling (DeepVariant 'min-base-quality' parameter).", true, 10);
$parser->addFloat("min_mq", "Minimum mapping quality used for variant calling (DeepVariant 'min-mapping-quality' parameter).", true, 20);
$parser->addFlag("test", "Run pipeline in test mode (e.g. use hard-coded cohort for methylation analysis)");
extract($parser->parse($argv));

// create logfile in output folder if no filepath is provided:
if (!file_exists($folder))
{
	exec2("mkdir -p $folder");
}

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($folder."/analyze_longread_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//determine processing system
$sys = load_system($system, $name);
$build = $sys['build'];
$platform = get_longread_sequencing_platform($sys['name_short']);

if ($test) trigger_error("Pipeline running in test mode!", E_USER_WARNING);

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

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
	if (in_array("re", $steps) && $rc_svs_exist)
	{
		trigger_error("Skipping step 're' - Report configuration with REs exists in NGSD!", E_USER_NOTICE);
		if (($key = array_search("re", $steps)) !== false) unset($steps[$key]);
	}
}

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$build);
}

$genome = genome_fasta($build);

//output file names:
//mapping
$bam_file = $folder."/".$name.".bam";
$cram_file = $folder."/".$name.".cram";
$used_bam_or_cram = ""; //BAM/CRAM file used for calling etc. This is a local tmp file if mapping was done and a file in the output folder if no mapping was done
// $unmapped_bam_file = $folder."/".$name.".mod.unmapped.bam";
$lowcov_file = $folder."/".$name."_".$sys["name_short"]."_lowcov.bed";
//methylation 
$modkit_track = $folder."/".$name."_modkit_track.bed.gz";
$modkit_summary = $folder."/".$name."_modkit_summary.txt";
$epigen_tsv = $folder."/".$name."_epigen.tsv";
//variant calling
$vcf_file = $folder."/".$name."_var.vcf.gz";
$vcf_file_annotated = $folder."/".$name."_var_annotated.vcf.gz";
$var_file = $folder."/".$name.".GSvar";
$roh_file = $folder."/".$name."_rohs.tsv";
$baf_file = $folder."/".$name."_bafs.igv";
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
//methylation files
$methyl_regions = repository_basedir()."/data/methylation/methylartist_catalog_grch38.tsv";
$methylation_table = $folder."/".$name."_var_methylation.tsv";
//db import
$qc_fastq  = $folder."/".$name."_stats_fastq.qcML";
$qc_map  = $folder."/".$name."_stats_map.qcML";
$qc_vc  = $folder."/".$name."_stats_vc.qcML";
$qc_other  = $folder."/".$name."_stats_other.qcML";

$name_sample_ps = explode("_", $name, 2);
$sample_name = $name_sample_ps[0];

if(!$sys['target_file']) trigger_error("Target region is missing in processing system config.", E_USER_NOTICE);

//check if target region covers whole genome (based on ROI size because pipeline test does not contain all chromosomes)
$check_chrs = bed_size(realpath($sys['target_file'])) > 3e9;
if(!$check_chrs) trigger_error("Target region does not cover whole genome. Cannot check for missing chromosomes in variant calls.", E_USER_NOTICE);

//mapping
if (in_array("ma", $steps))
{
	//determine input FASTQ and BAM files
	$unmapped_pattern = "{$sample_name}_??.mod.unmapped.bam";
	$unmapped_bam_files = glob("{$folder}/{$unmapped_pattern}");
	$old_bam_files = array_merge(glob("{$folder}/{$sample_name}_??.bam"), glob("{$folder}/{$sample_name}_??.cram"));
	$fastq_pattern = "{$sample_name}_??*.fastq.gz";
	$fastq_files = glob("{$folder}/{$fastq_pattern}");
	// preference:
	// 1. mod unmapped bam
	// 2. regular bam/cram
	// 3. fastq	
	if ((count($unmapped_bam_files) === 0) && (count($old_bam_files) === 0) && (count($fastq_files) === 0))
	{
		trigger_error("Found no input read files in BAM or FASTQ format!\nExpected file names are: unmapped BAM '{$unmapped_pattern}', mapped bam '{$sample_name}_??.[bam|cram]' or FASTQs '{$fastq_pattern}'", E_USER_ERROR);
	}
	
	// run mapping
	$mapping_minimap_options = [
		"-sample {$name}",
		"-threads {$threads}",
		"-system {$system}",
		"-qc_fastq {$qc_fastq}",
		"-qc_map {$qc_map}",
		"-softclip_supplements"	];
	if ($bam_output)
	{
		$mapping_minimap_options[] = "-out {$bam_file}";
		$mapping_minimap_options[] = "-bam_output";
		$used_bam_or_cram = $bam_file;
	}
	else
	{
		$mapping_minimap_options[] = "-out {$cram_file}";
		//create separate tmp folder with correct BamFile name to prevent problems with tools which uses this name as sample name:
		$local_bam = $parser->tempFolder("local_bam")."/".$name.".bam"; 
		$mapping_minimap_options[] = "-local_bam {$local_bam}";
		$used_bam_or_cram = $local_bam;
	}
	if (count($unmapped_bam_files) > 0)
	{
		$mapping_minimap_options[] = "-in_bam " . implode(" ", $unmapped_bam_files);	
	}
	elseif (count($old_bam_files) > 0)
	{
		//move BAMs/CRAMs to subfolder
		$old_bam_files_moved = array();
		$mapping_folder = "{$folder}/bams_for_mapping/";
		$parser->exec("mkdir", $mapping_folder);

		foreach ($old_bam_files as $old_bam) 
		{
			$basename = basename($old_bam);
			$parser->moveFile($old_bam, $mapping_folder.$basename);
			$old_bam_files_moved[] = $mapping_folder.$basename;
		}
		$mapping_minimap_options[] = "-in_bam " . implode(" ", $old_bam_files_moved);
	}
	else
	{
		$mapping_minimap_options[] = "-in_fastq " . implode(" ", $fastq_files);
	}
	
	$parser->execTool("Tools/mapping_minimap.php", implode(" ", $mapping_minimap_options));

	// create methylation track
	if (contains_methylation($used_bam_or_cram, 100, $build))
	{
		$args = array();
		$args[] = "-bam ".$used_bam_or_cram;
		$args[] = "-bed ".$modkit_track;
		$args[] = "-summary ".$modkit_summary;
		$args[] = "-threads ".$threads;
		$args[] = "-build ".$build;
		$parser->execTool("Tools/vc_modkit.php", implode(" ", $args));
	}

	//low-coverage report
	$parser->execApptainer("ngs-bits", "BedLowCoverage", "-in ".realpath($sys['target_file'])." -bam {$used_bam_or_cram} -out $lowcov_file -cutoff 20 -threads {$threads} -ref {$genome}", [$sys['target_file'], $genome, $used_bam_or_cram], [$folder]);
	if (db_is_enabled("NGSD"))
	{
		$parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", [$folder]);
	}

	//check if BAM/CRAM exist and have a appropriete size
	$bam_exists = file_exists($bam_file) && file_exists($bam_file.".bai"); 
	$cram_exists = file_exists($cram_file) && file_exists($cram_file.".crai"); 

	//get mapped read count
	if ($cram_exists) $rc_mapped_bam = get_read_count($cram_file, max(8, $threads), array("-F", "2304"), $sys['build']);
	else if ($bam_exists) $rc_mapped_bam = get_read_count($bam_file, max(8, $threads), array("-F", "2304"), $sys['build']);
	else trigger_error("Cannot delete FASTQ/BAM files - no BAM/CRAM file found!", E_USER_ERROR);
	trigger_error("Mapped BAM file read counts: ".$rc_mapped_bam, E_USER_NOTICE);

	//if reanalysis is performed: remove temporary subfolder
	if (isset($old_bam_files_moved) && count($old_bam_files_moved) > 0)
	{
		//check read counts of bams
		$rc_old = 0;
		foreach ($old_bam_files_moved as $old_bam) 
		{
			$read_count = get_read_count($old_bam, max(8, $threads), array("-F", "2304"), $sys['build']);
			trigger_error("Previous BAM/CRAM file read counts (".$old_bam."): ".$read_count, E_USER_NOTICE);
			$rc_old += $read_count;
		}

		//calculate relative difference
		$diff = abs($rc_old - $rc_mapped_bam);
		$rel_diff = $diff /(($rc_old + $rc_mapped_bam)/2);

		if (($rc_mapped_bam == $rc_old) || ($rel_diff < 0.0001))
		{
			// remove old BAM(s)
			if ($rc_mapped_bam == $rc_old) trigger_error("Read count of mapped and unmapped BAM(s) match. Deleting old BAM(s)/CRAM(s)...", E_USER_NOTICE);
			if ($rel_diff < 0.0001) trigger_error("Read count of new and old BAM(s)/CRAM(s) in allowed tolerance (<0.01%) (".($rel_diff*100)."%). Deleting old BAM(s)/CRAM(s)...", E_USER_NOTICE);
			$mapping_folder = "{$folder}/bams_for_mapping";
			$parser->exec("rm", "-r {$mapping_folder}");
			trigger_error("previous BAM/CRAM file removed!", E_USER_NOTICE);
		}
		else
		{
			trigger_error("Cannot delete old BAM/CRAM file(s) - mapped BAM file read counts doesn't match old BAM(s)/CRAM(s) (".($rel_diff*100)."%).", E_USER_ERROR);
		}
	}

	if (get_path("delete_fastq_files"))
	{
		//check if project overwrites the settings
		$preserve_fastqs = false;
		if (db_is_enabled("NGSD"))
		{
			$db = DB::getInstance("NGSD", false);
			$info = get_processed_sample_info($db, $name, false);
			if (!is_null($info))
			{
				$preserve_fastqs = $info['preserve_fastqs'];
			}
		}
		
		
		if(!$preserve_fastqs)
		{
			//get read count of input files
			if (count($unmapped_bam_files) > 0)
			{
				//check read counts of bams
				$rc_unmapped = 0;
				foreach ($unmapped_bam_files as $unmapped_bam_file) 
				{
					$read_count = get_read_count($unmapped_bam_file, max(8, $threads), array(), $sys['build']);
					trigger_error("Unmapped BAM file read counts (".$unmapped_bam_file."): ".$read_count, E_USER_NOTICE);
					$rc_unmapped += $read_count;
				}

				//calculate relative difference
				$diff = abs($rc_unmapped - $rc_mapped_bam);
				$rel_diff = $diff /(($rc_unmapped + $rc_mapped_bam)/2);

				if (($rc_mapped_bam == $rc_unmapped) || ($rel_diff < 0.0001))
				{
					// remove unmapped BAM(s)
					if ($rc_mapped_bam == $rc_unmapped) trigger_error("Read count of mapped and unmapped BAM(s) match. Deleting unmapped BAM(s)...", E_USER_NOTICE);
					if ($rel_diff < 0.0001) trigger_error("Read count of mapped and unmapped BAM(s) in allowed tolerance (<0.01%) (".($rel_diff*100)."%). Deleting unmapped BAM(s)...", E_USER_NOTICE);
					foreach ($unmapped_bam_files as $unmapped_bam_file) 
					{
						unlink($unmapped_bam_file);
						trigger_error("Unmapped BAM file removed!", E_USER_NOTICE);
					}
				}
				else
				{
					trigger_error("Cannot delete unmapped BAM file(s) - mapped BAM file read counts doesn't match unmapped BAM(s)(".($rel_diff*100)."%).", E_USER_ERROR);
				}

			}
			elseif (count($fastq_files) > 0)
			{
				//check read counts of bams
				$rc_unmapped = 0;
				foreach ($fastq_files as $fastq_file) 
				{
					$rc_pipeline = array();
					$rc_pipeline[] = array("zcat", $fastq_file);
					$rc_pipeline[] = array("wc", "-l");
					list($stdout, $stderr) = $parser->execPipeline($rc_pipeline, "FastQ read count", true);
					$read_count = ((int) $stdout[0]) / 4;
					trigger_error("FASTQ file read counts (".$fastq_file."): ".$read_count, E_USER_NOTICE);
					$rc_unmapped += $read_count;
				}

				//calculate relative difference
				$diff = abs($rc_unmapped - $rc_mapped_bam);
				$rel_diff = $diff /(($rc_unmapped + $rc_mapped_bam)/2);

				if (($rc_mapped_bam == $rc_unmapped) || ($rel_diff < 0.0001))
				{
					// remove unmapped BAM(s)
					if ($rc_mapped_bam == $rc_unmapped) trigger_error("Read count of mapped and FastQ(s) match. Deleting unmapped BAM(s)...", E_USER_NOTICE);
					if ($rel_diff < 0.0001) trigger_error("Read count of mapped BAM and FastQ(s) in allows tolerance (<0.01%) (".($rel_diff*100)."%). Deleting FastQ(s)...", E_USER_NOTICE);
					foreach ($fastq_files as $fastq_file) 
					{
						unlink($fastq_file);
						trigger_error("FastQ file removed!", E_USER_NOTICE);
					}
				}
				else
				{
					trigger_error("Cannot delete unmapped BAM file(s) - mapped BAM file read counts doesn't match unmapped BAM(s)(".($rel_diff*100)."%).", E_USER_ERROR);
				}
			}
			else
			{
				trigger_error("No FastQ/unmapped BAMs found!", E_USER_WARNING);
			}
				
		}
	}
}
else
{
	//set BAM/CRAM to use
	$used_bam_or_cram = file_exists($bam_file) ? $bam_file : $cram_file;
	$bam_output = file_exists($bam_file);

	if(!file_exists($used_bam_or_cram)) trigger_error("BAM/CRAM file not found!", E_USER_ERROR);
	
	//check genome build of BAM
	check_genome_build($used_bam_or_cram, $build);
}


//check gender after mapping
if(db_is_enabled("NGSD") && !$no_gender_check)
{
	$parser->execTool("Tools/db_check_gender.php", "-in $used_bam_or_cram -pid $name");	
}

//variant calling
if (in_array("vc", $steps))
{
	if ($platform == "PB")
	{
		$args = [];
		$args[] = "-model_type PACBIO";
		$args[] = "-bam ".$used_bam_or_cram;
		$args[] = "-out ".$vcf_file;
		$args[] = "-build ".$build;
		$args[] = "-threads ".$threads;
		$args[] = "-target ".$sys['target_file'];
		$args[] = "-target_extend 200";
		$args[] = "-min_af ".$min_af;
		$args[] = "-min_mq ".$min_mq;
		$args[] = "-min_bq ".$min_bq;
		$args[] = "-add_sample_header";
		$args[] = "-name ".$name;

		if ($gpu)
		{
			$args[] = "-gpu";
		}

		$parser->execTool("Tools/vc_deepvariant.php", implode(" ", $args));
	}
	else
	{

		//determine basecall model
		$basecall_model = get_basecall_model($used_bam_or_cram);
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
			else if ($sys["name_short"] == "LR-ONT-SQK-LSK114-SUP")
			{
				$basecall_model_path = get_path("clair3_models")."/r1041_e82_400bps_sup_v500/";
			}
			else
			{
				trigger_error("Unsupported processing system '".$sys["name_short"]."' provided!", E_USER_ERROR);
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
		$args[] = "-bam ".$used_bam_or_cram;
		$args[] = "-folder ".$folder;
		$args[] = "-name ".$name;
		$args[] = "-target ".$sys['target_file'];
		$args[] = "-target_extend 200";
		$args[] = "-threads ".$threads;
		$args[] = "-build ".$build;
		$args[] = "--log ".$parser->getLogFile();
		$args[] = "-model ".$basecall_model_path;
		if ($gpu) $args[] = "-gpu";
		
		$parser->execTool("Tools/vc_clair.php", implode(" ", $args));
	}

	//create b-allele frequency file
	$params = array();
	$params[] = "-vcf {$vcf_file}";
	$params[] = "-name {$name}";
	$params[] = "-out {$baf_file}";
	$params[] = "-downsample 100";
	$parser->execTool("Tools/baf_germline.php", implode(" ", $params));
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
				["", $parser->execApptainer("ngs-bits", "BedChunk", "-in ".realpath($sys['target_file'])." -n {$cnv_bin_size}", [$sys['target_file']], [], true)],
				["", $parser->execApptainer("ngs-bits", "BedAnnotateGC", "-clear -ref ".$genome, [$genome], [], true)],
				["", $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-out {$bed}", [], [dirname($bed)], true)]
			];
		$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
	}

	//create coverage profile
	$tmp_folder = $parser->tempFolder();
	$cov_file = $cov_folder."/{$name}.cov.gz";
	$cov_tmp_unzipped = $tmp_folder."/{$name}.cov";
	$cov_tmp = $cov_tmp_unzipped.".gz";
	$parser->execApptainer("ngs-bits", "BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam {$used_bam_or_cram} -in {$bed} -out {$cov_tmp_unzipped} -threads {$threads} -ref {$genome}", [$folder, $bed, $genome]);
	$parser->exec("gzip", "-9 {$cov_tmp_unzipped}");
	
	//copy coverage file to reference folder if valid
	if (db_is_enabled("NGSD") && is_valid_ref_sample_for_cnv_analysis($name))	{
		$parser->log("Moving coverage file to reference folder...");
		$parser->moveFile($cov_tmp, $cov_file);	}
	else
	{
		$cov_file = $cov_tmp;
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
	
	$parser->execTool("Tools/vc_clincnv_germline.php", implode(" ", $args), true);
	
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
	$parser->execTool("Tools/vc_sniffles.php", "-bam {$used_bam_or_cram} -include_mosaic -sample_ids {$name} -out {$sv_vcf_file} -threads {$threads} -build {$build} --log ".$parser->getLogFile());
}

//phasing
if (in_array("ph", $steps))
{
	//replace contigs in VCF header (incorrect sorted contigs lead to errors in CRAM file)
	$tmp_vcf = $parser->tempFile("_var.vcf");
	$contig_pipeline = array();
	$contig_pipeline[] = array("zcat", $vcf_file);
	$contig_pipeline[] = array("egrep", "-v \"^##contig=\" > {$tmp_vcf}");
	$parser->execPipeline($contig_pipeline, "contig removal");
	add_missing_contigs_to_vcf($sys['build'], $tmp_vcf);
	$parser->execApptainer("htslib", "bgzip", "-c {$tmp_vcf} > {$vcf_file}", [], [dirname($vcf_file)]);
	$parser->execApptainer("htslib", "tabix", "-f -p vcf {$vcf_file}", [], [dirname($vcf_file)]);
	
	//check for methylation
	$contains_methylation = contains_methylation($used_bam_or_cram);

	//create modcall file
	if($contains_methylation)
	{
		//delete prev modcall file
		if(file_exists($vcf_modcall)) $parser->exec("rm", $vcf_modcall);
		$args = array();
		$args[] = "modcall";
		$args[] = "-b {$used_bam_or_cram}";
		$args[] = "-r {$genome}";
		$args[] = "-t {$threads}";
		$args[] = "-o ".substr($vcf_modcall, 0, -4);
		$parser->execApptainer("longphase", "longphase", implode(" ", $args), [$folder,  $genome]);

		trigger_error("Methylation annotation detected. Using intermediate modcall step.", E_USER_NOTICE);
	}

	//run phasing by LongPhase on VCF files
	$phased_tmp = $parser->tempFile(".vcf", "longphase");
	$phased_sv_tmp = substr($phased_tmp,0,-4)."_SV.vcf";

	$args = array();
	$args[] = "phase";
	$args[] = "-s {$vcf_file}";
	$args[] = "-b {$used_bam_or_cram}";
	$args[] = "-r {$genome}";
	$args[] = "-t {$threads}";
	$args[] = "-o ".substr($phased_tmp, 0, -4);
	$args[] = $platform == "PB" ? "--pb" : "--ont";
	$args[] = "--indels";
	if (file_exists($sv_vcf_file))
	{
		$args[] = "--sv-file {$sv_vcf_file}";
	} 
	if ($contains_methylation) 
	{
		$args[] = "--mod-file {$vcf_modcall}";
	}
	$parser->execApptainer("longphase", "longphase", implode(" ", $args), [$folder, $genome], [$folder]);
	
	//create compressed file and index
	$parser->execApptainer("htslib", "bgzip", "-c $phased_tmp > {$vcf_file}", [], [dirname($vcf_file)]);
	$parser->execApptainer("htslib", "tabix", "-f -p vcf $vcf_file", [], [dirname($vcf_file)]);
	if (file_exists($sv_vcf_file))
	{
		$parser->execApptainer("htslib", "bgzip", "-c $phased_sv_tmp > {$sv_vcf_file}", [], [dirname($sv_vcf_file)]);
		$parser->execApptainer("htslib", "tabix", "-f -p vcf $sv_vcf_file", [], [dirname($sv_vcf_file)]);
	}

	//tag BAM file 
	$args = array();
	$args[] = "haplotag";
	$args[] = "-s {$vcf_file}";
	$args[] = "-b {$used_bam_or_cram}";
	$args[] = "-r {$genome}";
	$args[] = "-t {$threads}";
	if (file_exists($sv_vcf_file))
	{
		$args[] = "--sv-file {$sv_vcf_file}";
		$in_files[] = $sv_vcf_file;
	} 
	
	$tagged_bam_file = $parser->tempFile(".tagged.bam");
	$args[] = "-o ".dirname($tagged_bam_file)."/".basename2($tagged_bam_file);

	$parser->execApptainer("longphase", "longphase", implode(" ", $args), [$folder, $genome]);
	$parser->indexBam($tagged_bam_file, $threads);

	//check if read counts of tagged and untagged BAMs match
	compare_bam_read_count($used_bam_or_cram, $tagged_bam_file, max(8, $threads), true, true, 0.0, array(), $sys['build']);

	//use tagged bam in /tmp for further analysis
	$used_bam_or_cram = $tagged_bam_file;

	if ($bam_output)
	{
		//move BAM file in case copy fails:
		$parser->moveFile($bam_file, $bam_file.".backup.bam");
		$parser->moveFile($bam_file.".bai", $bam_file.".backup.bam.bai");

		//use BAM
		$parser->copyFile($tagged_bam_file, $bam_file);
		$parser->copyFile($tagged_bam_file.".bai", $bam_file.".bai");

		// check if copy was successful
		if (!file_exists($bam_file) || filesize($tagged_bam_file) != filesize($bam_file))
		{
			trigger_error("Error during coping BAM file! File sizes don't match!", E_USER_ERROR);
		}

		//delete backup file
		unlink($bam_file.".backup.bam");
		unlink($bam_file.".backup.bam.bai");
	}
	else
	{
		//move BAM file in case conversion fails:
		$parser->moveFile($cram_file, $cram_file.".backup.cram");
		$parser->moveFile($cram_file.".crai", $cram_file.".backup.cram.crai");

		//convert bam to cram
		$parser->execTool("Tools/bam_to_cram.php", "-bam {$tagged_bam_file} -cram {$cram_file} -build ".$sys['build']." -threads {$threads}");
		$parser->indexBam($cram_file, $threads);
		
		//check if read counts of tagged and untagged CRAMs match
		compare_bam_read_count($cram_file.".backup.cram", $cram_file, max(8, $threads), true, true, 0.0, array(), $sys['build']);
		
		//delete backup file
		unlink($cram_file.".backup.cram");
		unlink($cram_file.".backup.cram.crai");
	}
}
else if(in_array("ma", $steps) || in_array("vc", $steps) || in_array("sv", $steps))
{
	trigger_error("Mapping or (SV) variant calling done without phasing step! Output files might not be phased correctly!", E_USER_WARNING);
}

// repeat expansion
if (in_array("re", $steps))
{
	//Repeat-expansion calling using straglr
	$variant_catalog = repository_basedir()."/data/repeat_expansions/straglr_variant_catalog_grch38.bed";
	$parser->execTool("Tools/vc_straglr.php", "-include_partials -in {$used_bam_or_cram} -out {$straglr_file} -loci {$variant_catalog} -threads {$threads} -build {$build} --log ".$parser->getLogFile());
}


// paralogous gene calling
if (in_array("pg", $steps))
{
	$parser->execTool("Tools/vc_paraphase.php", "-folder {$folder} -name {$name} ".(ends_with($used_bam_or_cram, ".bam")?"-local_bam {$used_bam_or_cram} ":"")."-build {$build} -threads {$threads} --log ".$parser->getLogFile());
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
		$parser->execTool("Tools/annotate.php", implode(" ", $args));

		//check for truncated output
		if ($check_chrs) $parser->execTool("Tools/check_for_missing_chromosomes.php", "-in {$vcf_file_annotated} -max_missing_perc 5");
		
		//ROH detection
		$in_files = [];
		$in_files[] = $folder;
		$in_files[] = repository_basedir()."/data/gene_lists/genes.bed";
		$in_files[] = repository_basedir()."/data/misc/roh_exclude_regions.bed";
		$args = [];
		$args[] = "-in $vcf_file_annotated";
		$args[] = "-out $roh_file";
		$args[] = "-var_min_dp 15";
		$args[] = "-var_min_q 15";
		$args[] = "-var_af_keys gnomADg_AF";
		$args[] = "-exclude ".repository_basedir()."/data/misc/roh_exclude_regions.bed";
		$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file)) $in_files[] = $omim_file;
		$args[] = "-annotate ".repository_basedir()."/data/gene_lists/genes.bed ".(file_exists($omim_file) ? $omim_file : "");
		$parser->execApptainer("ngs-bits", "RohHunter", implode(" ", $args), $in_files);

		//PRS calculation 
		$prs_folder = repository_basedir()."/data/misc/prs/";
		$prs_scoring_files = glob($prs_folder."/*_".$build.".vcf");
		if (count($prs_scoring_files) > 0)
		{
			$parser->execApptainer("ngs-bits", "VcfCalculatePRS", "-in {$vcf_file} -bam {$used_bam_or_cram} -out $prs_file -prs ".implode(" ", $prs_scoring_files)." -ref $genome -long_read", [$folder, $genome, $prs_folder]);
		}

		//determine ancestry
		if (ngsbits_build($build) != "non_human")
		{
			$parser->execApptainer("ngs-bits", "SampleAncestry", "-in {$vcf_file} -out {$ancestry_file} -build ".ngsbits_build($build), [$folder]);
		}
	}

	// annotate CNV file
	if (file_exists($cnv_file))
	{
		$repository_basedir = repository_basedir();
		$data_folder = get_path("data_folder");
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$cnv_file}", [$folder, "{$repository_basedir}/data/misc/af_genomes_imgag.bed"]);
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$cnv_file}", [$folder, "{$repository_basedir}/data/misc/cn_pathogenic.bed"]);
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed -no_duplicates -url_decode -out {$cnv_file}", [$folder, "{$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed"]);
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2025-09.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$cnv_file}", [$folder, "{$data_folder}/dbs/ClinVar/clinvar_cnvs_2025-09.bed"]);

		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2025_2.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$cnv_file}", [$folder, $hgmd_file]);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$omim_file} -no_duplicates -url_decode -out {$cnv_file}", [$folder, $omim_file]);
		}

		//annotate gene info
		if(db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "CnvGeneAnnotation", "-in {$cnv_file} -add_simple_gene_names -out {$cnv_file}", [$folder]);
		}
		
		//annotate overlap with pathogenic CNVs
		if (db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "NGSDAnnotateCNV", "-in {$cnv_file} -out {$cnv_file}", [$folder]);
		}

		//check for truncated output
		if ($check_chrs) $parser->execTool("Tools/check_for_missing_chromosomes.php", "-in {$cnv_file}");
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
		$parser->execApptainer("ngs-bits", "VcfToBedpe", "-in $sv_vcf_file -out $bedpe_file", [$folder]);

		//add gene info annotation
		if (db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "BedpeGeneAnnotation", "-in $bedpe_file -out $bedpe_file -add_simple_gene_names", [$folder]);;
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
			$args = array(
				"-in $bedpe_file",
				"-out $bedpe_file",
				"-processing_system ".$sys["name_short"],
				"-ann_folder $ngsd_annotation_folder",
			);

			if (db_is_enabled("NGSD"))
			{
				$db = DB::getInstance("NGSD", false);
				$ps_id = get_processed_sample_id($db, $name, false);
	
				if ($ps_id != -1)
				{
					$args[] = "-ps_name $name";
				}
				else 
				{
					trigger_error("No processed sample ID found for sample ".$name.", skipping count annotation by disease group!", E_USER_WARNING);
				}
			}

			$parser->execApptainer("ngs-bits", "BedpeAnnotateCounts", implode(" ", $args), [$folder, $ngsd_annotation_folder]);

			$sys_specific_density_file = $ngsd_annotation_folder."sv_breakpoint_density_".$sys["name_short"].".igv";
			if (file_exists($sys_specific_density_file))
			{
				$parser->execApptainer("ngs-bits", "BedpeAnnotateBreakpointDensity", "-in {$bedpe_file} -out {$bedpe_file} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv -density_sys {$sys_specific_density_file}", [$folder, $ngsd_annotation_folder]);
			}
			else
			{
				$parser->execApptainer("ngs-bits", "BedpeAnnotateBreakpointDensity", "-in {$bedpe_file} -out {$bedpe_file} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv", [$folder, $ngsd_annotation_folder]);
			}
			

			// check if files changed during annotation
			foreach ($ngsd_sv_files as $filename)
			{
				if ($db_file_dates[$filename] != filemtime($ngsd_annotation_folder.$filename))
				{
					trigger_error("Annotation file '".$ngsd_annotation_folder.$filename."' has changed during annotation!",E_USER_ERROR);
				}
			}

			//annotate class 4 and 5 pathogenic SVs
			if (db_is_enabled("NGSD"))
			{
				$parser->execApptainer("ngs-bits", "NGSDAnnotateSV", "-in {$bedpe_file} -out {$bedpe_file}", [$folder]);
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
			$parser->execApptainer("ngs-bits", "BedpeAnnotateFromBed", "-in $bedpe_file -out $bedpe_file -bed $omim_file -url_decode -replace_underscore -col_name OMIM", [$folder, $omim_file]);
		}

		//add CNV overlap annotation
		if (file_exists($cnv_file))
		{
			$parser->execApptainer("ngs-bits", "BedpeAnnotateCnvOverlap", "-in $bedpe_file -out $bedpe_file -cnv $cnv_file", [$folder]);
		}

		//write genotype in own column
		$parser->execApptainer("ngs-bits", "BedpeExtractGenotype", "-in $bedpe_file -out $bedpe_file -include_unphased", [$folder]);

		//extract columns
		$parser->execApptainer("ngs-bits", "BedpeExtractInfoField", "-in $bedpe_file -out $bedpe_file -info_fields SVLEN,SUPPORT,COVERAGE,AF", [$folder]);

		//update sample entry 
		update_gsvar_sample_header($bedpe_file, array($name=>"Affected"));	
		
		//check for truncated output
		if ($check_chrs) $parser->execTool("Tools/check_for_missing_chromosomes.php", "-in {$bedpe_file}");
	}
}

// INFO: done after annotation to have phasing track available
// methylation calling
if (in_array("me", $steps))
{
	if (!contains_methylation($used_bam_or_cram)) trigger_error("BAM file doesn't contain methylation info! Skipping step 'me'", E_USER_WARNING);
	else 
	{
		$args = [];
		$args[] = "-folder {$folder}";
		$args[] = "-name {$name}";
		if (ends_with($used_bam_or_cram, ".bam")) $args[] = "-local_bam {$used_bam_or_cram}";
		$args[] = "-out {$methylation_table}";
		$args[] = "-build {$build}";
		$args[] = "-threads {$threads}";
		$args[] = "--log ".$parser->getLogFile();
		if ($test)
		{
			$args[] = "-regions ".repository_basedir()."/test/data/create_methyl_plot_regions_pipeline_test.tsv";
			$args[] = "-custom_cohort_table ".repository_basedir()."/test/data/create_methyl_plot_custom_cohort.tsv";
			$args[] = "-test";
		}
		else
		{
			$args[] = "-regions {$methyl_regions}";
		}
		$parser->execTool("Tools/create_methyl_plot.php", implode(" ", $args));
		
		// create Epigen TSV file
		$epic_id_file = get_path("data_folder")."/dbs/illumina-epicids/EPIC-8v2-0_A1.csv";
		$parser->execApptainer("ngs-bits", "BedToEpigen", "-in {$modkit_track} -out {$epigen_tsv} -id_file {$epic_id_file} -sample {$name}", [$epic_id_file, $modkit_track], [dirname($epigen_tsv)]);
		
	}
}

// collect other QC terms - if CNV or SV calling was done
if ((in_array("ma", $steps)) || (in_array("cn", $steps) || in_array("sv", $steps) || in_array("an", $steps)))
{
	$terms = [];
	$sources = [];

	//Basecall model
	if (file_exists($used_bam_or_cram))
	{
		$basecall_model = strtr(get_basecall_model($used_bam_or_cram), "_", " ");
		if ($basecall_model != "")
		{
			$terms[] = "QC:2000149\t{$basecall_model}";
			$sources[] = $used_bam_or_cram;
		}	
	}
	
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
				if ($ll>=12) # changed to 12 for long-reads
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
		$parser->execApptainer("ngs-bits", "TsvToQC", "-in $tmp -out $qc_other -sources ".implode(" ", $sources), [$folder]);
	}
}

// create Circos plot - if small variant, CNV or SV calling was done
if (in_array("an", $steps))
{
	if (file_exists($cnv_file))
	{
		if (file_exists($cnv_file2))
		{
			$parser->execTool("Tools/create_circos_plot.php", "-folder $folder -name $name -lrgs -build ".$build);
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
	if(file_exists($ancestry_file)) $parser->execTool("Tools/db_import_ancestry.php", "-id {$name} -in {$ancestry_file}");
	
	//import QC
	$qc_files = array($qc_fastq, $qc_map);
	if (file_exists($qc_vc)) $qc_files[] = $qc_vc; 
	if (file_exists($qc_other)) $qc_files[] = $qc_other; 
	$parser->execApptainer("ngs-bits", "NGSDImportSampleQC", "-ps $name -files ".implode(" ", $qc_files)." -force", [$folder]);

	if (!$skip_wgs_check)
	{
		//check similarity to (sr) WGS/WES sample
		list($stdout, $stderr, $exit_code) = $parser->execApptainer("ngs-bits", "NGSDSameSample", "-ps $name -system_type WGS,WES");
		//parse stdout
		$same_samples = array();
		foreach ($stdout as $line) 
		{
			if (starts_with($line, "#")) continue;
			if (trim($line) == "") continue;
			$parts = explode("\t", $line);
			$same_samples[] = trim($parts[0]);
		}
		if (count($same_samples) < 1)
		{
			trigger_error("No related WGS/WES sample found for {$name}. Cannot perform similarity check!", E_USER_WARNING);
		} 
		else
		{
			//check same-sample correlation for each WGS sample
			$min_corr = 0.85;

			foreach ($same_samples as $processed_sample) 
			{
				$ps_gsvar = trim($parser->execApptainer("ngs-bits", "SamplePath", "-ps {$processed_sample} -type GSVAR")[0][0]);
				if (!file_exists($ps_gsvar))
				{
					trigger_error("GSvar file {$ps_gsvar} not found! Skipping sample similarity check", E_USER_WARNING);
					continue;
				}
				$output = $parser->execApptainer("ngs-bits", "SampleSimilarity", "-in {$var_file} {$ps_gsvar} -mode gsvar -build ".ngsbits_build($sys['build']), [$folder, $ps_gsvar]);
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
	$args = [];
	if (file_exists($var_file))
	{
		//check genome build
		check_genome_build($var_file, $build);
		
		$args[] = "-var {$var_file}";
	}
	if (file_exists($cnv_file))
	{
		//check genome build
		//this is not possible for CNVs because the file does not contain any information about it
		
		$args[] = "-cnv {$cnv_file}";
	}
	if (file_exists($bedpe_file))
	{
		//check genome build
		check_genome_build($bedpe_file, $build);
		
		$args[] = "-sv {$bedpe_file}";
	}
	if (file_exists($straglr_file))
	{
		//check genome build
		check_genome_build($straglr_file, $build);
		
		$args[] = "-re {$straglr_file}";
	}
	if (count($args)>0)
	{
		$args[] = "-ps {$name}";
		$parser->execApptainer("ngs-bits", "NGSDAddVariantsGermline", implode(" ", $args), [$folder]);
	}
}

?>
