<?php

/**
	@page analyze_cfdna
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("analyze_cfdna", "cfDNA analysis pipeline.");
$parser->addInfile("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);
//optional
$parser->addString("tumor_id", "Related tumor processed sample.", true, "");
$parser->addInfile("tumor_bam", "BAM file of related tumor processed sample.", true, "");
$parser->addInfile("target", "Custom target region as BED file.", true, "");
$parser->addInfile("monitoring_vcf", "VCF containing patient-specific variants (and IDs).", true);
$parser->addFlag("skip_tumor", "Skip comparison with related tumor sample");
$parser->addInt("base_extend", "Number of bases the target region is extended (default: 60)", true, 60);
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'name' is a valid processed sample name).", true);
$steps_all = array("ma", "vc", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, db=import into NGSD.", true, "ma,vc,db");
$parser->addFlag("annotation_only", "Performs only a re-annotation of the already created variant calls.");
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFloat("min_corr", "The minimum sample genotype correlation which is used for tumor-cfDNA comparison (default: 0.80)", true, 0.80);
$parser->addInt("min_mapq", "The minimum mapping quality for reads to be considered.", true, 0);
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
$parser->addFlag("no_post_filter", "Use the unfiltered umiVar output to generate the GSvar file.");
extract($parser->parse($argv));

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($folder."/analyze_cfdna_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//TODO determine values
//low coverage cutoff
$lowcov_cutoff = 100;

//log parameters
$parser->log("Parameter:");
$parser->log("   folder          = {$folder}");
$parser->log("   name            = {$name}");
$parser->log("   tumor_id        = {$tumor_id}");
$parser->log("   tumor_bam       = {$tumor_bam}");
$parser->log("   target          = {$target}");
$parser->log("   monitoring_vcf  = {$monitoring_vcf}");
$parser->log("   skip_tumor      = {$skip_tumor}");
$parser->log("   base_extend     = {$base_extend}");
$parser->log("   system          = {$system}");
$parser->log("   steps           = {$steps}");
$parser->log("   annotation_only = {$annotation_only}");
$parser->log("   threads         = {$threads}");
$parser->log("   min_corr        = {$min_corr}");

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//determine processing system
$sys = load_system($system, $name);
$genome = genome_fasta($sys['build']);

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);
}

// determine analysis type
if ($sys['type']=="cfDNA (patient-specific)" || $sys['type']=="cfDNA")
{
	$is_patient_specific = $sys['type']=="cfDNA (patient-specific)";
}
else
{
	trigger_error("Unsupported system type '".$sys['type']."'! Pipeline only supports cfDNA samples", E_USER_ERROR);	
}

// for std cfDNA -> use processing system bed as default target
if (!$is_patient_specific && $target == "")
{
	// get target region from processing system
	$target = realpath($sys['target_file']);
}

//database
if (db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
}

// overwrite tumor_bam and tumor_id if run with -skip_tumor
if ($skip_tumor)
{
	$tumor_bam = "";
	$tumor_id = "";
}


if (isset($db))
{
	// get tumor id from NGSD
	if (!$skip_tumor && ($tumor_id == ""))
	{
		//related samples
		list($sample_name, $ps_num) = explode("_", $name);
		$res_samples = $db->executeQuery("SELECT id, name FROM sample WHERE name=:name", ["name" => $sample_name]);
		if (count($res_samples) !== 1)
		{
			trigger_error("Could not find sample for processed sample {$name}!", E_USER_WARNING);
		}
		$sample_id = $res_samples[0]['id'];
		$res = $db->executeQuery("SELECT * FROM sample_relations WHERE relation='tumor-cfDNA' AND (sample1_id=:sid OR sample2_id=:sid)", array("sid" => $sample_id));

		$psamples = [];
		$ps_ids = [];
		foreach ($res as $row)
		{
			$sample_id_annotation = $row['sample1_id'] != $sample_id ? $row['sample1_id'] : $row['sample2_id'];
			$res = $db->executeQuery("SELECT ps.sample_id, ps.process_id, ps.processing_system_id, ps.quality, sys.id, sys.type, CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as psample, ps.id FROM processed_sample as ps, processing_system as sys, sample as s WHERE ps.sample_id=:sid AND ps.processing_system_id=sys.id AND ps.sample_id=s.id AND sys.type!='cfDNA (patient-specific)' ORDER BY ps.process_id ASC", array("sid" => $sample_id_annotation));
			$psamples = array_merge($psamples, array_column($res, 'psample'));
		}
		if (count($psamples) > 1)
		{
			trigger_error("Found more than one possible referenced tumor: " . implode(", ", $psamples), E_USER_WARNING);

			//identify tumor sample by cfDNA panel
			$sys_id = $db->getValue("SELECT id FROM processing_system WHERE name_short='".$sys["name_short"]."'");

			$panel_counts = 0;
			foreach ($psamples as $tsample) 
			{
				$tumor_ps_id = get_processed_sample_id($db, $tsample);
				$res = $db->executeQuery("SELECT id FROM cfdna_panels WHERE vcf IS NOT NULL AND vcf <> '' AND processing_system_id=:sysid AND tumor_id=:tid", array("sysid" => $sys_id, "tid" => $tumor_ps_id));
				$panel_counts = count($res);
				if ($panel_counts > 0)
				{
					$tumor_id = $tsample;
					trigger_error("Tumor id with cfDNA panel extracted from the NGSD: ${tumor_id}", E_USER_NOTICE);
					break;
				}	
			}

			// fallback if no cfDNA panel was found
			if ($tumor_id == "")
			{
				$tumor_id = $psamples[0];
				trigger_error("Found more than one referenced tumor, using first one: " . implode(" ", $psamples), E_USER_WARNING);
			}
			
		}
		else if (count($psamples) === 0)
		{
			trigger_error("Could not find any related tumor processed sample!", E_USER_WARNING);
		}
		else
		{
			$tumor_id = $psamples[0];
			trigger_error("Tumor id extracted from the NGSD: ${tumor_id}", E_USER_NOTICE);
		}
	}

	// get tumor bam file from NGSD
	if (!$annotation_only && !$skip_tumor && ($tumor_id != "") && ($tumor_bam == ""))
	{
		$psinfo_tumor = get_processed_sample_info($db, $tumor_id);
		$tumor_bam = $psinfo_tumor['ps_bam'];
	}


	// get related cfDNA samples from NGSD
	if (!$skip_tumor && ($tumor_id != ""))
	{

		list($tumor_name, ) = explode("_", $tumor_id);
		$res_samples = $db->executeQuery("SELECT id, name FROM sample WHERE name=:name", ["name" => $tumor_name]);
		if (count($res_samples) !== 1)
		{
			trigger_error("Could not find sample for selected tumor sample {$tumor_id}!", E_USER_WARNING);
		}
		$tumor_db_id = $res_samples[0]['id'];
		$res = $db->executeQuery("SELECT * FROM sample_relations WHERE relation='tumor-cfDNA' AND sample1_id=:sid", array("sid" => $tumor_db_id));

		$related_cfdna_samples = array();
		foreach ($res as $row)
		{
			$res = $db->executeQuery("SELECT ps.sample_id, ps.process_id, ps.processing_system_id, ps.quality, sys.id, sys.type, CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as psample, ps.id FROM processed_sample as ps, processing_system as sys, sample as s WHERE ps.sample_id=:sid AND ps.processing_system_id=sys.id AND ps.sample_id=s.id AND (sys.type='cfDNA (patient-specific)' OR sys.type='cfDNA') ORDER BY ps.process_id ASC", array("sid" => $row['sample2_id']));
			$related_cfdna_samples = array_merge($related_cfdna_samples, array_column($res, 'psample'));
		}

		// remove current sample from list of related samples
		if (($key = array_search($name, $related_cfdna_samples)) !== false) 
		{
			unset($related_cfdna_samples[$key]);
		}

		// get BAMs and GSvar files (for post-filtering)
		$related_cfdna_bams = array();
		$related_cfdna_gsvars = array();
		foreach ($related_cfdna_samples as $sample) 
		{
			$sample_info = get_processed_sample_info($db, $sample);
			
			//GSvar files
			$gsvar = $sample_info['ps_folder'].$sample_info['ps_name'].".GSvar";
			if (file_exists($gsvar))
			{
				$related_cfdna_gsvars[$sample] = $gsvar;
			}
			else
			{
				trigger_error("GSvar file for sample $sample not found!", E_USER_WARNING);
			}

			//BAM files
			$bam = $sample_info['ps_bam'];
			if (file_exists($bam))
			{
				$related_cfdna_bams[$sample] = $bam;
			}
			else
			{
				trigger_error("BAM file for sample $sample not found!", E_USER_WARNING);
			}
		}

		trigger_error("Related cfDNA samples: ".implode(", ", $related_cfdna_samples), E_USER_NOTICE);

	}

	// generate bed files for calling (only required for mapping and VC)
	if (!$annotation_only && (in_array("ma", $steps) || in_array("vc", $steps)))
	{
		if ($is_patient_specific)
		{
			if (($target == "") || ($monitoring_vcf == ""))
			{
				// get patient-specific target region from NGSD
				if ($tumor_id == "") trigger_error("Valid tumor id is required to get patient specific target region!", E_USER_ERROR);

				$tumor_ps_id = get_processed_sample_id($db, $tumor_id);
				$sys_id = $db->getValue("SELECT id FROM processing_system WHERE name_short='".$sys["name_short"]."'");

				if ($target == "")
				{
					$target = "{$folder}/{$name}_target_region.bed";
					$bed_content = $db->getValue("SELECT bed FROM cfdna_panels WHERE tumor_id=${tumor_ps_id} AND processing_system_id=${sys_id}");
					file_put_contents($target, $bed_content);
				}
				else
				{
					copy($target, "{$folder}/{$name}_target_region.bed");
				}
				
				if ($monitoring_vcf == "")
				{
					$monitoring_vcf = "{$folder}/{$name}_cfdna_panel.vcf";
					$vcf_content = $db->getValue("SELECT vcf FROM cfdna_panels WHERE tumor_id=${tumor_ps_id} AND processing_system_id=${sys_id}");
					file_put_contents($monitoring_vcf, $vcf_content);
				}
				else
				{
					copy($monitoring_vcf, "{$folder}/{$name}_cfdna_panel.vcf");
				}
			}
		}
	}
}
else
{
	trigger_error("WARNING: No connection to the NGSD! All additional information that is stored in the database will be missing in the output!", E_USER_WARNING);
}

//output file names:
//mapping
$bamfile = "{$folder}/{$name}.bam";
$bamfile_raw = "{$folder}/{$name}_before_dedup.bam";
$lowcov_file = "{$folder}/{$name}_{$sys['name_short']}-{$tumor_id}_lowcov.bed";
//variant calling
$vcffile = "{$folder}/{$name}_var.vcf";
//db import
$qc_fastq = "{$folder}/{$name}_stats_fastq.qcML";
$qc_map  = "{$folder}/{$name}_stats_map.qcML";
$qc_cfdna = "{$folder}/{$name}_stats_cfDNA.qcML";

//alternative mrd calculation
$bg_mrd = "{$folder}/umiVar/{$name}_bg.mrd";
$bg_mrd_unfiltered = "{$folder}/umiVar/{$name}_bg_unfiltered.mrd";
$bg_monitoring = "{$folder}/umiVar/{$name}_monitoring_counts.tsv";
$bg_monitoring_unfiltered = "{$folder}/umiVar/{$name}_monitoring_counts_unfiltered.tsv";

//check if all required info are available
if (!$annotation_only && (in_array("ma", $steps) || in_array("vc", $steps)))
{
	if ($is_patient_specific && !$skip_tumor)
	{
		if ($tumor_id == "") trigger_error("Required parameter 'tumor_id' missing! Can not perform analysis.", E_USER_ERROR);
		if ($tumor_bam == "") trigger_error("Required parameter 'tumor_bam' missing! Can not perform analysis.", E_USER_ERROR);
	}

	if ($target == "") trigger_error("No target region defined! Can not perform analysis.", E_USER_ERROR);
	if ($is_patient_specific && ($monitoring_vcf == "")) trigger_error("No monitoring VCF defined! Can not perform analysis.", E_USER_ERROR);
}

// for annotation_only: check if all files are available
if ($annotation_only)
{
	if(in_array("vc", $steps) && !file_exists($vcffile))
	{
		trigger_error("VCF for reannotation is missing. Skipping 'vc' step!", E_USER_WARNING);
		if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
	}  
}


//extend target region for calling
if (!$annotation_only && (in_array("ma", $steps) || in_array("vc", $steps)))
{
	$target_extended = $parser->tempFile("_extended.bed");
	$pipeline = [];
	if ($base_extend > 0)
	{
		$pipeline[] = ["", $parser->execApptainer("ngs-bits", "BedExtend", "-in {$target} -n {$base_extend}", [$target], [], true)];
	}
	else
	{
		$pipeline[] = [ "cat", $target];
	}
	$pipeline[] = ["", $parser->execApptainer("ngs-bits", "BedSort", "", [], [], true)];
	$pipeline[] = ["", $parser->execApptainer("ngs-bits", "BedMerge", "-out {$target_extended}", [], [], true)];
	$parser->execPipeline($pipeline, "extend target region");
}

// for annotation_only: check if all files are available
if ($annotation_only)
{
	if(in_array("vc", $steps) && !file_exists($vcffile))
	{
		trigger_error("VCF for reannotation is missing. Skipping 'vc' step!", E_USER_WARNING);
		if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
	} 
}

//mapping
if (in_array("ma", $steps))
{
	//determine input FASTQ files
	$in_for = "{$folder}/*_R1_00?.fastq.gz";
	$in_rev = "{$folder}/*_R2_00?.fastq.gz";
	$in_index = "{$folder}/*_index_*.fastq.gz";
	
	$files1 = glob($in_for);
	$files2 = glob($in_rev);
	$files_index = glob($in_index);
	if (count($files1) != count($files2))
	{
		trigger_error("Found mismatching forward and reverse read file count!\n Forward: $in_for\n Reverse: $in_rev.", E_USER_ERROR);
	}
	if (!empty($files_index) && count($files_index) != count($files1))
	{
		trigger_error("Found mismatching index read file count!", E_USER_ERROR);
	}
	if (empty($files1))
	{
		trigger_error("No FastQ files found in sample folder!", E_USER_ERROR);
	}
	
	$args = [
		"-in_for ".implode(" ", $files1),
		"-in_rev ".implode(" ", $files2),
		"-clip_overlap",
		"-system " .$system,
		"-out_folder ".$folder,
		"-out_name ".$name,
		"-threads ".$threads,
		"--log ".$parser->getLogFile(),
		"-bam_output" //TODO Leon: CRAM does not work in pipeline test - error message is "[E::process_one_read] CIGAR and query sequence are of different length"
	];
	if ($min_mapq > 0 ) $args[] = "-min_mapq $min_mapq";
	if (!empty($files_index)) $args[] = "-in_index " . implode(" ", $files_index);
	$parser->execTool("Tools/mapping.php", implode(" ", $args));

	//low-coverage report, based on patient specific positions
	$parser->execApptainer("ngs-bits", "BedLowCoverage", "-in $target -bam $bamfile -out $lowcov_file -cutoff $lowcov_cutoff -threads {$threads} -ref {$genome}", [$target, $folder, $genome]);
	if (db_is_enabled("NGSD"))
	{
		$parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", [$folder]);
	}

	$parser->execApptainer("ngs-bits", "MappingQC", "-cfdna -roi {$target_extended} -in {$bamfile} -out {$qc_map} -ref {$genome} -build ".ngsbits_build($sys['build']), [$folder, $genome]);

}

//variant calling
if (in_array("vc", $steps))
{
	// skip VC if only annotation should be done
	if (!$annotation_only)
	{
		// check if reference genome fits BAM file
		check_genome_build($bamfile, $sys['build']);

		$args = [
			"-bam", $bamfile,
			"-target", $target_extended,
			"-build", $sys['build'],
			"-folder", $folder."/umiVar",
			"--log", $parser->getLogFile()
		];

		// //TODO: model still needed?
		// $model = get_path("data_folder")."/dbs/cfdna_caller/{$sys['name_short']}.txt";
		// if (file_exists($model))
		// {
		// 	$args[] = "-model {$model}";
		// }

		if ($is_patient_specific)
		{
			$args[] = "-monitoring_vcf ${monitoring_vcf}";
		}
		$parser->execTool("Tools/vc_cfdna.php", implode(" ", $args));

		//copy vcf (monitoring variants for monitoring / filtered for normal cfDNA)
		if ($is_patient_specific)
		{
			$parser->copyFile("{$folder}/umiVar/{$name}_monitoring.vcf", $vcffile);
		}
		else
		{
			if ($no_post_filter)
			{
				$parser->copyFile("{$folder}/umiVar/{$name}.vcf", $vcffile);
			}
			else
			{
				$parser->copyFile("{$folder}/umiVar/{$name}_hq.vcf", $vcffile);
			}
			
		}

		// sort VCF
		$parser->execApptainer("ngs-bits", "VcfSort", "-in $vcffile -out $vcffile", [$folder]);
		
		// mark off-target reads
		$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in $vcffile -mark off-target -reg $target -out $vcffile", [$folder, $target]);

		// validate VCF
		$parser->execApptainer("ngs-bits", "VcfCheck", "-in $vcffile -lines 0 -ref {$genome}", [$folder, $genome]);


		// cfDNA QC
		$in_files = array();
		$args = [
			"-bam", $bamfile,
			"-cfdna_panel", $target,
			"-out", $qc_cfdna,
			"-build", ngsbits_build($sys['build']),
			"-ref", $genome
		];
		$in_files[] = $folder;
		$in_files[] = $target;
		$in_files[] = $genome;
		// add tumor BAM
		if ($tumor_bam != "")
		{
			$args[] = "-tumor_bam ".$tumor_bam;
			$in_files[] = $tumor_bam;
		}
		// add related BAMs
		if (isset($related_cfdna_bams) && count($related_cfdna_bams) > 0)
		{
			$args[] = "-related_bams ".implode(" ", array_values($related_cfdna_bams));
			foreach (array_values($related_cfdna_bams) as $in_file) 
			{
				$in_files[] = $in_file;
			}
		}
		// add umiVar error rates
		$error_rate_file = "{$folder}/umiVar/error_rates.txt";
		if (file_exists($error_rate_file))
		{
			$args[] = "-error_rates ".$error_rate_file;
		}

		//run CfDnaQC
		$parser->execApptainer("ngs-bits", "CfDnaQC", implode(" ", $args), $in_files);
	}

	// check if reference genome fits VCF file
	if ($annotation_only) check_genome_build($vcffile, $sys['build']);

	// annotate VCF 
	$vcffile_annotated = "$folder/${name}_var_annotated.vcf.gz";
	$parser->execTool("Tools/annotate.php", "-out_name $name -out_folder $folder -vcf $vcffile -somatic -threads $threads -system $system");

	//add sample info to VCF header
	$tmp_vcf = $parser->tempFile(".vcf");
	$s = Matrix::fromTSV($vcffile_annotated);
	$comments = $s->getComments();
	$comments[] = gsvar_sample_header($name, array("IsTumor" => "yes"), "#", "");
	$s->setComments(sort_vcf_comments($comments));
	$s->toTSV($tmp_vcf);

	// zip and index vcf file
	$parser->execApptainer("htslib", "bgzip", "-c $tmp_vcf > $vcffile_annotated", [], [dirname($vcffile_annotated)]);
	$parser->execApptainer("htslib", "tabix", "-f -p vcf $vcffile_annotated", [], [dirname($vcffile_annotated)]);

	// convert vcf to GSvar
	$gsvar_file = "$folder/${name}.GSvar";
	$args = array("-in $tmp_vcf", "-out $gsvar_file", "-t_col $name", "-cfdna");
	$parser->execTool("Tools/vcf2gsvar_somatic.php", implode(" ", $args));

	if (!$skip_tumor && ($tumor_id != ""))
	{
		//add normal AF & counts
		if (isset($db))
		{
			// get normal sample
			$psinfo_tumor = get_processed_sample_info($db, $tumor_id);
			if ($psinfo_tumor['normal_name'] == "") 
			{
				trigger_error("Normal id not found!", E_USER_WARNING);
			}
			else
			{
				$psinfo_normal = get_processed_sample_info($db, $psinfo_tumor['normal_name']);
				// annotate AF and depth
				$in_files = array();
				$in_files[] = $psinfo_normal['ps_bam'];
				$out_files = array();
				$out_files[] = $folder;
				$parser->execApptainer("ngs-bits", "VariantAnnotateFrequency", "-in $gsvar_file -out $gsvar_file -depth -name normal -bam ".$psinfo_normal['ps_bam'], $in_files, $out_files);
			}
		}
		else
		{
			trigger_error("No NGSD connection! Can not annotate normal AF!", E_USER_WARNING);
		}
	}

	// perform GSvar post-filtering
	if ($is_patient_specific)
	{
		if (!isset($related_cfdna_gsvars)) $related_cfdna_gsvars = array();
		$related_cfdna_gsvars[$name] = $gsvar_file;
		ksort($related_cfdna_gsvars); //sort by sample name --> //TODO: sort by sample id

		// modify output file list so that only the GSvar file of the current sample is written
		$related_cfdna_gsvars_output = array_fill(0, count($related_cfdna_gsvars), "");
		$sample_pos = array_search($name, array_keys($related_cfdna_gsvars));
		$related_cfdna_gsvars_output[$sample_pos] = $gsvar_file;
		
		// initialize input and output file list for umivar container execution
		$in_files = array();
		$in_files[] = $gsvar_file;

		$temp_logfile = $parser->tempFile(".log");
		$args = [
			implode(",", $related_cfdna_gsvars),
			implode(",", $related_cfdna_gsvars_output),
			"--log_file ".$temp_logfile
		];

		$in_files = array_merge($in_files, $related_cfdna_gsvars);
		
		//get tumor-normal GSvar
		if ((isset($db)) && ($tumor_id != ""))
		{
			$tumor_gsvar_files_db = $db->getValues("SELECT `gsvar_file` FROM `secondary_analysis` WHERE `type` = 'somatic' AND `gsvar_file` LIKE '%{$tumor_id}%' ORDER BY `type` ASC ");
			$tumor_gsvar_files = array();
			foreach($tumor_gsvar_files_db as $tn_gsvar_file)
			{
				if (file_exists($tn_gsvar_file))
				{
					$tumor_gsvar_files[] = $tn_gsvar_file;
				}
				else
				{
					trigger_error("Tumor-normal GSvar file $tn_gsvar_file not found!", E_USER_WARNING);
				}
			}
			
			if (count($tumor_gsvar_files) > 0)
			{
				$args[] = "--tumor_samples ".implode(",", $tumor_gsvar_files);
				$in_files = array_merge($in_files, $tumor_gsvar_files);
			}
			else
			{
				trigger_error("No tumor-normal analysis for sample {$tumor_id} found!", E_USER_WARNING);
			}
		}

		//post-filtering
		$parser->execApptainer("umiVar", "cfDNA_postfiltering.py", implode(" ", $args), $in_files, [$folder]);
		$parser->log("post-filtering log: ", file($temp_logfile));

	}

	// calculate alternative MRD & monitoring counts
	if ($is_patient_specific)
	{

		//filtered output
		$in_files = array();
		$in_files[] = $folder;
		$in_files[] = $monitoring_vcf;
		$args = array();
		$args[] = $folder."/umiVar";
		$args[] = $monitoring_vcf;
		$args[] = $bg_mrd;
		$args[] = $bg_monitoring;
		$parser->execApptainer("umiVar", "calculateMRD.py", implode(" ", $args), $in_files);

		//unfiltered output
		$args = array();
		$args[] = $folder."/umiVar";
		$args[] = $monitoring_vcf;
		$args[] = $bg_mrd_unfiltered;
		$args[] = $bg_monitoring_unfiltered;
		$args[] = "--max_af 1.5";
		$args[] = "--keep_gonosomes";
		$args[] = "--keep_indels";
		$parser->execApptainer("umiVar", "calculateMRD.py", implode(" ", $args), $in_files);
	}
	
}

//check sample similarity with referenced tumor
if (!($annotation_only || $skip_tumor))
{
	if (($tumor_bam != "") && (in_array("ma", $steps) || in_array("vc", $steps)))
	{
		$output = $parser->execApptainer("ngs-bits", "SampleSimilarity", "-in {$bamfile} {$tumor_bam} -mode bam -ref {$genome} -build ".ngsbits_build($sys['build']), [$folder, $tumor_bam, $genome]);
		$correlation = explode("\t", $output[0][1])[3];
		if ($correlation < $min_corr)
		{
			trigger_error("The genotype correlation of cfDNA and tumor ({$tumor_id}) is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
		}
		else
		{
			trigger_error("The genotype correlation of cfDNA and tumor ({$tumor_id}) is {$correlation}.", E_USER_NOTICE);
		}

		// calculate similarity between related cfDNA samples
		foreach ($related_cfdna_bams as $cfdna_sample => $cfdna_bam) 
		{
			$output = $parser->execApptainer("ngs-bits", "SampleSimilarity", "-in {$bamfile} {$cfdna_bam} -mode bam -ref {$genome} -build ".ngsbits_build($sys['build']), [$folder, $cfdna_bam, $genome]);
			$correlation = explode("\t", $output[0][1])[3];
			trigger_error("The genotype correlation of cfDNA and related sample ({$cfdna_sample}) is {$correlation}.", E_USER_NOTICE);
		}
	}
	else
	{
		trigger_error("Skipping similarity check!", E_USER_WARNING);
	}
}

//import to database
if (in_array("db", $steps))
{
	$qc_files = array($qc_fastq, $qc_map);
	if (file_exists($qc_cfdna))
	{
		$qc_files[] = $qc_cfdna;
	}
	//import QC
	$parser->execApptainer("ngs-bits", "NGSDImportSampleQC", "-ps {$name} -files ".implode(" ", $qc_files)." -force", [$folder]);

	//check gender
	$parser->execTool("Tools/db_check_gender.php", "-in {$bamfile} -pid {$name} --log ".$parser->getLogFile());	
}

?>
