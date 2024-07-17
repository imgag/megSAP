<?php

/**
	@page analyze
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("analyze", "Complete NGS analysis pipeline.");
$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);
//optional
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'name' is a valid processed sample name).", true);
$steps_all = array("ma", "vc", "cn", "sv", "re", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, cn=copy-number analysis, sv=structural-variant analysis, db=import into NGSD.", true, "ma,vc,cn,sv,re,db");
$parser->addFloat("min_af", "Minimum VAF cutoff used for variant calling (freebayes 'min-alternate-fraction' parameter).", true, 0.1);
$parser->addFloat("min_bq", "Minimum base quality used for variant calling (freebayes 'min-base-quality' parameter).", true, 10);
$parser->addFloat("min_mq", "Minimum mapping quality used for variant calling (freebayes 'min-mapping-quality' parameter).", true, 20);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("clip_overlap", "Soft-clip overlapping read pairs.");
$parser->addFlag("no_abra", "Skip realignment with ABRA.");
$parser->addFlag("no_trim", "Skip adapter trimming with SeqPurge.");
$parser->addFlag("correction_n", "Use Ns for errors by barcode correction.");
$parser->addFlag("somatic", "Set somatic single sample analysis options (i.e. correction_n, clip_overlap).");
$parser->addFlag("annotation_only", "Performs only a reannotation of the already created variant calls.");
$parser->addFlag("use_dragen", "Use Illumina DRAGEN server for mapping, small variant and structural variant calling.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
$parser->addString("rna_sample", "Processed sample name of the RNA sample which should be used for annotation.", true, "");
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
$is_wes = $sys['type']=="WES";
$is_wgs = $sys['type']=="WGS";
$is_panel = $sys['type']=="Panel" || $sys['type']=="Panel Haloplex";
$is_wgs_shallow = $sys['type']=="WGS (shallow)";
$has_roi = $sys['target_file']!="";
$build = $sys['build'];
$genome = genome_fasta($build);

//handle somatic flag
if ($somatic)
{
	$clip_overlap = true;
	$correction_n = true;
}

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

//checks in case DRAGEN should be used
if ($use_dragen)
{
	if (!in_array("ma", $steps) && in_array("vc", $steps) && !file_exists($folder."/dragen_variant_calls/{$name}_dragen.vcf.gz")) 
	{
		trigger_error("DRAGEN small variant calls have to be present in the folder {$folder}/dragen_variant_calls for the use of DRAGEN without the mapping step!", E_USER_ERROR);
	}
	if (!in_array("ma", $steps) && in_array("sv", $steps) && get_path("use_dragen_sv_calling") && !file_exists($folder."/dragen_variant_calls/{$name}_dragen_svs.vcf.gz")) 
	{
		trigger_error("DRAGEN structural variant calls have to be present in the folder {$folder}/dragen_variant_calls for the use of DRAGEN without the mapping step!", E_USER_ERROR);
	}
}

//remove invalid steps
if (in_array("cn", $steps) && !$has_roi)
{
	trigger_error("Skipping step 'cn' - Copy number analysis is only supported for processing systems with target region BED file!", E_USER_NOTICE);
	if (($key = array_search("cn", $steps)) !== false) unset($steps[$key]);
}
if (in_array("sv", $steps) && $is_wgs_shallow)
{
	trigger_error("Skipping step 'sv' - Structural variant calling is not supported for shallow WGS samples!", E_USER_NOTICE);
	if (($key = array_search("sv", $steps)) !== false) unset($steps[$key]);
}

if (db_is_enabled("NGSD") && !$annotation_only)
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
$bamfile = $folder."/".$name.".bam";
$cramfile = $folder."/".$name.".cram";
$used_bam_or_cram = ""; //BAM/CRAM file used for calling etc. This is a local tmp file if mapping was done and a file in the output folder if no mapping was done
$lowcov_file = $folder."/".$name."_".$sys["name_short"]."_lowcov.bed";
$somatic_custom_panel = get_path("data_folder") . "/enrichment/somatic_VirtualPanel_v5.bed";
//variant calling
$vcffile = $folder."/".$name."_var.vcf.gz";
$vcffile_annotated = $folder."/".$name."_var_annotated.vcf.gz";
$varfile = $folder."/".$name.".GSvar";
$rohfile = $folder."/".$name."_rohs.tsv";
$baffile = $folder."/".$name."_bafs.igv";
$ancestry_file = $folder."/".$name."_ancestry.tsv";
$prsfile = $folder."/".$name."_prs.tsv";
//copy-number calling
$cnvfile = $folder."/".$name."_cnvs_clincnv.tsv";
$cnvfile2 = $folder."/".$name."_cnvs_clincnv.seg";
//structural variant calling
$sv_manta_file = $folder ."/". $name . "_manta_var_structural.vcf.gz";
$bedpe_out = substr($sv_manta_file,0,-6)."bedpe";
//repeat expansions
$expansion_hunter_file = $folder."/".$name."_repeats_expansionhunter.vcf";
//db import
$qc_fastq  = $folder."/".$name."_stats_fastq.qcML";
$qc_map  = $folder."/".$name."_stats_map.qcML";
$qc_vc  = $folder."/".$name."_stats_vc.qcML";
$qc_other  = $folder."/".$name."_stats_other.qcML";

// for annotation_only: check if all files are available
if ($annotation_only)
{
	if(in_array("vc", $steps) && !file_exists($vcffile))
	{
		trigger_error("VCF for reannotation is missing. Skipping 'vc' step!", E_USER_WARNING);
		if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
	} 

	if(in_array("cn", $steps) && !file_exists($cnvfile))
	{
		trigger_error("CN file for reannotation is missing. Skipping 'cn' step!", E_USER_WARNING);
		if (($key = array_search("cn", $steps)) !== false) unset($steps[$key]);
	} 

	if(in_array("sv", $steps) && !file_exists($bedpe_out))
	{
		trigger_error("SV file for reannotation is missing. Skipping 'sv' step!", E_USER_WARNING);
		if (($key = array_search("sv", $steps)) !== false) unset($steps[$key]);
	} 
}

//mapping
if (in_array("ma", $steps))
{
	//find FASTQ input files
	$in_for = $folder."/*_R1_00?.fastq.gz";
	$files1 = glob($in_for);
	$in_rev = $folder."/*_R2_00?.fastq.gz";
	$files2 = glob($in_rev);
	$in_index = $folder."/*_index_*.fastq.gz";
	$files_index = glob($in_index);
	
	//check number of FASTQs for forward/reverse is equal
	if (count($files1)!=count($files2))
	{
		trigger_error("Found mismatching forward and reverse read file count!\n Forward: $in_for\n Reverse: $in_rev.", E_USER_ERROR);
	}
	
	//generate FASTQs from BAM/CRAM is missing
	if (count($files1)==0)
	{
		$source_file = "";
		$cram_genome = $genome;
		if(file_exists($bamfile))
		{
			$source_file = $bamfile;
		}
		else if (file_exists($cramfile))
		{
			$source_file = $cramfile;
			if (check_genome_build($cramfile, "GRCh38", false)!=1) //if genome build is not GRCh38 the CRAM has to be from SolveRD > use 1000G reference genome with decoy sequences and without 'chr' 
			{
				$cram_genome = get_path("data_folder")."/genomes/GRCh37.nochr.fa";
				trigger_error("Reference genome of CRAM is not GRCh38. Using SolveRD hg19 reference genome to convert CRAM to FASTQ: {$cram_genome}!", E_USER_NOTICE);
			}
		}
		if ($source_file=="")
		{
			trigger_error("Found no read files matching '$in_for' or '$in_rev' and no BAM ($bamfile) or CRAM file ($cramfile). Cannot perform mapping without raw data!", E_USER_ERROR);
		}
		
		trigger_error("No FASTQ files found in sample folder. Generating FASTQs from {$source_file}!", E_USER_NOTICE);

		// extract reads from BAM file
		$in_fq_for = $folder."/{$name}_BamToFastq_R1_001.fastq.gz";
		$in_fq_rev = $folder."/{$name}_BamToFastq_R2_001.fastq.gz";
		$tmp1 = $parser->tempFile(".fastq.gz");
		$tmp2 = $parser->tempFile(".fastq.gz");
		$parser->exec("{$ngsbits}BamToFastq", "-in $source_file -ref $cram_genome -out1 $tmp1 -out2 $tmp2", true);
		$parser->moveFile($tmp1, $in_fq_for);
		$parser->moveFile($tmp2, $in_fq_rev);
		
		// use generated fastq files for mapping
		$files1 = array($in_fq_for);
		$files2 = array($in_fq_rev);
	}
	
	$args = [];
	if($clip_overlap) $args[] = "-clip_overlap";
	if($no_abra) $args[] = "-no_abra";
	if($no_trim) $args[] = "-no_trim";
	if($correction_n) $args[] = "-correction_n";
	if(!empty($files_index)) $args[] = "-in_index " . implode(" ", $files_index);
	if($use_dragen) $args[] = "-use_dragen";
	if($somatic) $args[] = "-somatic_custom_map";
	$used_bam_or_cram = $parser->tempFolder("local_bam")."/".$name.".bam"; //local copy of BAM file to reduce IO over network when mapping is done 
	$parser->execTool("Pipelines/mapping.php", "-in_for ".implode(" ", $files1)." -in_rev ".implode(" ", $files2)." -system $system -out_folder $folder -out_name $name -local_bam $used_bam_or_cram ".implode(" ", $args)." -threads $threads");
	
	//low-coverage report
	if ($has_roi && !$is_wgs_shallow)
	{	
		$parser->exec("{$ngsbits}BedLowCoverage", "-in ".$sys['target_file']." -bam $used_bam_or_cram -out $lowcov_file -cutoff 20 -threads {$threads} -ref {$genome}", true);
		if (db_is_enabled("NGSD"))
		{
			$parser->exec("{$ngsbits}BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", true);
		}
	}

	//delete fastq files after mapping - remove files only if sample doesn't contain UMIs and the corresponding setting is set in the settings.ini
	if ($sys['umi_type']=="n/a" && get_path("delete_fastq_files")) 
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
			$fastq_files = array_merge($files1, $files2);
			if(count($fastq_files)>0)
			{
				//determine FASTQ size
				$fastq_file_size = 0;
				foreach($fastq_files as $fq_file)
				{
					$fastq_file_size += filesize($fq_file);
				}
			
				//check if BAM/CRAM exist and have a appropriete size
				$bam_exists = file_exists($bamfile) && file_exists($bamfile.".bai"); 
				$cram_exists = file_exists($cramfile) && file_exists($cramfile.".crai"); 
				if ($bam_exists)
				{
					if (filesize($bamfile) > 0.3 * $fastq_file_size)
					{
						foreach($fastq_files as $fq_file)
						{
							unlink($fq_file);
						}
					}
					else
					{
						trigger_error("Cannot delete FASTQ files - BAM file smaller than 30% of FASTQ.", E_USER_ERROR);
					}
				}
				else if ($cram_exists)
				{
					if (filesize($cramfile) > 0.1 * $fastq_file_size)
					{
						foreach($fastq_files as $fq_file)
						{
							unlink($fq_file);
						}
					}
					else
					{
						trigger_error("Cannot delete FASTQ files - CRAM file smaller than 10% of FASTQ.", E_USER_ERROR);
					}
				}
				else
				{
					trigger_error("Cannot delete FASTQ files - no BAM/CRAM file found!", E_USER_ERROR);
				}
			}
		}
	}
}
else if (file_exists($bamfile) || file_exists($cramfile))
{	
	//set BAM/CRAM to use
	$used_bam_or_cram = file_exists($bamfile) ? $bamfile : $cramfile;
	
	//check genome build of BAM
	check_genome_build($used_bam_or_cram, $build);

	//QC for already mapped/called samples from the NovaSeq X
	if(!file_exists($qc_map))
	{
		//QC
		$params = array("-in $used_bam_or_cram", "-out {$qc_map}", "-ref ".genome_fasta($sys['build']), "-build ".ngsbits_build($sys['build']));
		if ($sys['target_file']=="" || $sys['type']=="WGS" || $sys['type']=="WGS (shallow)")
		{
			$params[] = "-wgs";
		}
		else
		{
			$params[] = "-roi ".$sys['target_file'];
		}
		if ($sys['build']!="GRCh38")
		{
			$params[] = "-no_cont";
		}
		if ($somatic && file_exists($somatic_custom_panel))
		{
			$params[] = "-somatic_custom_bed $somatic_custom_panel";
		}
		if (!file_exists($qc_fastq))
		{
			$params[] = "-read_qc $qc_fastq";
		}
		$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $params), true);
	}	
	
	//low-coverage regions for already mapped/called samples from the NovaSeq X
	if(!file_exists($lowcov_file))
	{
		if ($has_roi && !$is_wgs_shallow)
		{	
			$parser->exec("{$ngsbits}BedLowCoverage", "-in ".$sys['target_file']." -bam $used_bam_or_cram -out $lowcov_file -cutoff 20 -threads {$threads}", true);
			if (db_is_enabled("NGSD"))
			{
				$parser->exec("{$ngsbits}BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", true);
			}
		}
	}
}

//variant calling
if (in_array("vc", $steps))
{
	// skip VC if only annotation should be done
	if (!$annotation_only)
	{		
		//Do not call standard pipeline if there is only mitochondiral chrMT in target region
		$only_mito_in_target_region = false;
		if ($has_roi) 
		{
			$only_mito_in_target_region = exec2("cat ".$sys['target_file']." | cut -f1 | uniq")[0][0] == "chrMT";
		}
		//activate mito-calling for sWGS
		if ($is_wgs_shallow) $only_mito_in_target_region = true;

		//perform main variant calling on autosomes/genosomes
		if(!$only_mito_in_target_region)
		{
			if ($use_dragen)
			{
				if (!in_array("ma", $steps)) trigger_error("'-use_dragen' with no mapping step provided. Using old DRAGEN VCF for small variant calling.", E_USER_NOTICE);
				$pipeline = [];
				$dragen_output_vcf = $folder."/dragen_variant_calls/{$name}_dragen.vcf.gz";
				$pipeline[] = array("zcat", $dragen_output_vcf);
				
				//filter by target region (extended by 200) and quality 5
				$target = $parser->tempFile("_roi_extended.bed");
				$parser->exec($ngsbits."BedExtend"," -in ".$sys['target_file']." -n 200 -out $target -fai ".$genome.".fai", true);
				$pipeline[] = array($ngsbits."VcfFilter", "-reg {$target} -qual 5 -filter_clear");

				//split multi-allelic variants
				$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");

				//normalize all variants and align INDELs to the left
				$pipeline[] = array($ngsbits."VcfLeftNormalize", "-stream -ref $genome");

				//sort variants by genomic position
				$pipeline[] = array($ngsbits."VcfStreamSort", "");

				//zip
				$pipeline[] = array("bgzip", "-c > $vcffile", false);

				//execute pipeline
				$parser->execPipeline($pipeline, "Dragen small variants post processing");

				//mark off-target variants
				$tmp = $parser->tempFile("_offtarget.vcf");
				$parser->exec($ngsbits."VariantFilterRegions", "-in $vcffile -mark off-target -reg ".$sys['target_file']." -out $tmp", true);
				
				//remove variants with less than 3 alternative observations
				$tmp2 = $parser->tempFile("_ad.vcf");
				$hr = gzopen2($tmp, "r");
				$hw = fopen2($tmp2, "w");
				while(!gzeof($hr))
				{
					$line = trim(gzgets($hr));
					if (strlen($line)==0) continue;
					if ($line[0]=="#")
					{
						fwrite($hw, $line."\n");
						continue;
					}
					
					//filter by AO
					$parts = explode("\t", $line); //chr, pos, id, ref, alt, qual, filter, info, format, sample
					$format = explode(":", $parts[8]);
					$sample = explode(":", $parts[9]);
					$sample = array_combine($format, $sample);
					if (isset($sample['DP']) && isset($sample['AF']))
					{
						if($sample['DP'] * $sample['AF']<2.9) continue; //not 3 because of rounding errors
					}
					
					fwrite($hw, $line."\n");
				}
				fclose($hr);
				fclose($hw);
				
				//bgzip
				$parser->exec("bgzip", "-c $tmp2 > $vcffile", false);

				//index output file
				$parser->exec("tabix", "-p vcf $vcffile", false); //no output logging, because Toolbase::extractVersion() does not return
			}
			else
			{
				$args = [];
				$args[] = "-bam ".$used_bam_or_cram;
				$args[] = "-out ".$vcffile;
				$args[] = "-build ".$build;
				$args[] = "-threads ".$threads;
				if ($has_roi)
				{
					$args[] = "-target ".$sys['target_file'];
					$args[] = "-target_extend 200";
				}
				$args[] = "-min_af ".$min_af;
				$args[] = "-min_mq ".$min_mq;
				$args[] = "-min_bq ".$min_bq;
				$parser->execTool("NGS/vc_freebayes.php", implode(" ", $args));
			}
		}
		
		//perform special variant calling for mitochondria
		$mito = enable_special_mito_vc($sys) || $only_mito_in_target_region;
		if ($mito)
		{
			$vcffile_mito = $parser->tempFile("_mito.vcf.gz");
			if ($use_dragen)
			{
				$pipeline = [];
				$dragen_output_vcf = $folder."/dragen_variant_calls/{$name}_dragen.vcf.gz";
				$pipeline[] = array("zcat", $dragen_output_vcf);
				
				//filter by target region and quality 5
				$pipeline[] = array($ngsbits."VcfFilter", "-reg chrMT:1-16569 -qual 5 -filter_clear");

				//split multi-allelic variants
				$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");

				//normalize all variants and align INDELs to the left
				$pipeline[] = array($ngsbits."VcfLeftNormalize", "-stream -ref $genome");

				//sort variants by genomic position
				$pipeline[] = array($ngsbits."VcfStreamSort", "");

				//zip
				$pipeline[] = array("bgzip", "-c > $vcffile_mito", false);

				//execute pipeline
				$parser->execPipeline($pipeline, "Dragen small variants post processing (chrMT)");

				//index output file
				$parser->exec("tabix", "-p vcf $vcffile_mito", false); //no output logging, because Toolbase::extractVersion() does not return
			}
			else
			{
				$target_mito = $parser->tempFile("_mito.bed");
				file_put_contents($target_mito, "chrMT\t0\t16569");
				
				$args = [];
				$args[] = "-in ".$used_bam_or_cram;
				$args[] = "-out ".$vcffile_mito;
				$args[] = "-build ".$build;
				$args[] = "-target ".$target_mito;
				$args[] = "-min_af 0.01";
				$args[] = "-max_af 1.00";
				$args[] = "-max_gnomad_af 1.00";
				$args[] = "-min_obs 2";
				$args[] = "-min_mq 20";
				$parser->execTool("NGS/vc_mosaic.php", implode(" ", $args));
			}
		
			if($only_mito_in_target_region) 
			{
				$parser->copyFile($vcffile_mito, $vcffile);
			}
		}
		
		//Add header to VCF file
		$hr = gzopen2($vcffile, "r");
		$vcf = $parser->tempFile("_unzipped.vcf");
		$hw = fopen2($vcf, "w");
		while(!gzeof($hr))
		{
			$line = trim(gzgets($hr));
			if (strlen($line)==0) continue;
			if ($line[0]=="#" && $line[1]!="#")
			{
				fwrite($hw, "##ANALYSISTYPE=GERMLINE_SINGLESAMPLE\n");
				fwrite($hw, "##PIPELINE=".repository_revision(true)."\n");
				fwrite($hw, gsvar_sample_header($name, array("DiseaseStatus"=>"affected")));
			}
			fwrite($hw, $line."\n");
		}
		fclose($hr);
		
		//Add mitochondrial variants to vcffile in case mito was called and it is not a pure mitochondrial sample
		if($mito && !$only_mito_in_target_region)
		{
			$hr = gzopen2($vcffile_mito, "r");
			while(!gzeof($hr))
			{
				$line = trim(gzgets($hr));
				if ($line=="" || $line[0]=="#") continue;
				fwrite($hw, $line."\n");
			}
			fclose($hr);
		}
		fclose($hw);
		
		//call low mappability variants
		if ($is_wgs || ($is_wes && $has_roi) || ($is_panel && $has_roi))
		{
			//determine region
			$mapq0_regions = repository_basedir()."data/misc/low_mappability_region/mapq_eq0.bed";
			if ($is_wgs)
			{
				$roi_low_mappabilty = $mapq0_regions;
			}
			else
			{
				$roi_low_mappabilty = $parser->tempFile("_mapq0.bed");
				$parser->exec("{$ngsbits}BedIntersect", "-in ".$sys['target_file']." -in2 $mapq0_regions -out $roi_low_mappabilty", true);
			}
			
			if (bed_size($roi_low_mappabilty)>0)
			{
				$tmp_low_mappability = $parser->tempFile("_low_mappability.vcf.gz");
				$args = [];
				$args[] = "-bam ".$used_bam_or_cram;
				$args[] = "-out ".$tmp_low_mappability;
				$args[] = "-build ".$build;
				$args[] = "-threads ".$threads;
				$args[] = "-target ".$roi_low_mappabilty;
				$args[] = "-min_af ".$min_af;
				$args[] = "-min_mq 0";
				$args[] = "-min_bq ".$min_bq;
				$parser->execTool("NGS/vc_freebayes.php", implode(" ", $args));
			
				//unzip
				$tmp_low_mappability2 = $parser->tempFile("_low_mappability.vcf");
				$parser->exec("zcat", "$tmp_low_mappability > $tmp_low_mappability2", true);
				
				//add to main variant list
				$tmp2 = $parser->tempFile("_merged_low_mappability.vcf");
				$parser->exec("{$ngsbits}VcfAdd", "-in $vcf -in2 $tmp_low_mappability2 -skip_duplicates -filter low_mappability -filter_desc Variants_in_reads_with_low_mapping_score. -out $tmp2", true);
				$parser->moveFile($tmp2, $vcf);
			}
		}
			
		//call mosaic variants on target region (WES or Panel) or exonic/splicing region (WGS)
		if (ngsbits_build($build)!="non_human" && ($is_wgs || ($is_wes && $has_roi) || ($is_panel && $has_roi)))
		{
			$tmp_mosaic = $parser->tempFile("_mosaic.vcf");
			$args = [];
			$args[] = "-in {$used_bam_or_cram}";
			$args[] = "-out {$tmp_mosaic}";
			$args[] = "-no_zip";
			$args[] = "-target ".(($is_panel || $is_wes) ? $sys['target_file'] : repository_basedir()."data/gene_lists/gene_exons_pad20.bed");
			$args[] = "-threads ".$threads;
			$args[] = "-build ".$build;
			$args[] = "-min_af 0.03";
			$args[] = "-min_obs ".($is_wgs ?  "1" : "2");	
			$parser->execTool("NGS/vc_mosaic.php", implode(" ", $args));
			
			//add to main variant list
			$tmp2 = $parser->tempFile("_merged_mosaic.vcf");
			$parser->exec("{$ngsbits}VcfAdd", "-in $vcf -in2 $tmp_mosaic -skip_duplicates -filter mosaic -filter_desc Putative_mosaic_variants. -out $tmp2", true);
			$parser->moveFile($tmp2, $vcf);
		}
				
		//sort and remove unused contig lines
		$parser->exec("{$ngsbits}VcfSort", "-in $vcf -remove_unused_contigs -out $vcf", true);
		
		//zip and index
		$parser->exec("bgzip", "-c $vcf > $vcffile", false); //no output logging, because Toolbase::extractVersion() does not return
		$parser->exec("tabix", "-f -p vcf $vcffile", false); //no output logging, because Toolbase::extractVersion() does not return
		
		//create b-allele frequency file
		$params = array();
		$params[] = "-vcf {$vcffile}";
		$params[] = "-name {$name}";
		$params[] = "-out {$baffile}";
		if ($is_wgs)
		{
			$params[] = "-downsample 100";
		}
		$parser->execTool("NGS/baf_germline.php", implode(" ", $params));
	}
	else
	{
		check_genome_build($vcffile, $build);
	}

	//annotation
	$args = [];
	$args[] = "-out_name ".$name;
	$args[] = "-out_folder ".$folder;
	$args[] = "-system ".$system;
	$args[] = "--log ".$parser->getLogFile();
	$args[] = "-threads ".$threads;
	if($rna_sample != "") $args[] = "-rna_sample ".$rna_sample;
	$parser->execTool("Pipelines/annotate.php", implode(" ", $args));

	//ROH detection
	if ($is_wes || $is_wgs)
	{
		$args = [];
		$args[] = "-in $vcffile_annotated";
		$args[] = "-out $rohfile";
		$args[] = "-var_af_keys gnomADg_AF";
		$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; //optional because of license
		$args[] = "-annotate ".repository_basedir()."/data/gene_lists/genes.bed ".(file_exists($omim_file) ? $omim_file : "");
		$parser->exec("{$ngsbits}RohHunter", implode(" ", $args), true);
	}

	//PRS calculation 
	if ($is_wgs)
	{
		$prs_folder = repository_basedir()."/data/misc/prs/";
		$prs_scoring_files = glob($prs_folder."/*_".$build.".vcf");
		if (count($prs_scoring_files) > 0)
		{
			$parser->exec("{$ngsbits}VcfCalculatePRS", "-in $vcffile -bam $used_bam_or_cram -out $prsfile -prs ".implode(" ", $prs_scoring_files)." -ref $genome", true);
		}
	}
	
	//determine ancestry
	if (ngsbits_build($build) != "non_human")
	{
		$parser->exec($ngsbits."SampleAncestry", "-in {$vcffile} -out {$ancestry_file} -build ".ngsbits_build($build), true);
	}
}

//copy-number analysis
if (in_array("cn", $steps))
{
	// skip CN calling if only annotation should be done
	if (!$annotation_only)
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
		//WGS: create folder for binned coverage data - if missing
		if ($is_wgs || $is_wgs_shallow)
		{

			$bin_size = get_path("cnv_bin_size_wgs");
			if ($is_wgs_shallow) $bin_size = get_path("cnv_bin_size_shallow_wgs");
			$bin_folder = "{$ref_folder}/bins{$bin_size}/";
			if (!is_dir($bin_folder))
			{

				mkdir($bin_folder);
				if (!chmod($bin_folder, 0777))
				{
					trigger_error("Could not change privileges of folder '{$bin_folder}'!", E_USER_ERROR);
				}
			}
			$cov_folder = $bin_folder;
		}

		
		//create BED file with GC and gene annotations - if missing
		if ($is_wgs || $is_wgs_shallow)
		{
			$bed = $ref_folder."/bins{$bin_size}.bed";
			if (!file_exists($bed))
			{
				$pipeline = [
						["{$ngsbits}BedChunk", "-in ".$sys['target_file']." -n {$bin_size}"],
						["{$ngsbits}BedAnnotateGC", "-clear -ref ".$genome],
						["{$ngsbits}BedAnnotateGenes", "-out {$bed}"]
					];
				$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
			}
		}
		else
		{
			$bed = $ref_folder."/roi_annotated.bed";
			if (!file_exists($bed))
			{
				$pipeline = [
						["{$ngsbits}BedAnnotateGC", "-in ".$sys['target_file']." -clear -ref ".$genome],
						["{$ngsbits}BedAnnotateGenes", "-out {$bed}"],
					];
				$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
			}
		}

		//create coverage profile
		$tmp_folder = $parser->tempFolder();
		$cov_file = $cov_folder."/{$name}.cov";
		if (!file_exists($cov_file) || filemtime($cov_file)<filemtime($used_bam_or_cram))
		{
			$parser->log("Calculating coverage file for CN calling...");
			$cov_tmp = $tmp_folder."/{$name}.cov";
			$parser->exec("{$ngsbits}BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam {$used_bam_or_cram} -in {$bed} -out {$cov_tmp} -threads {$threads} -ref {$genome}", true);
			
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
		);
		if ($is_wgs)
		{
			$args[] = "-max_cnvs 2000";
		}
		else if ($is_wgs_shallow)
		{
			$args[] = "-max_cnvs 400";
			$args[] = "-skip_super_recall";
			$args[] = "-regions 3";
		}
		else
		{
			$args[] = "-max_cnvs 200";
		}
		if($is_wgs || $is_wes || $is_wgs_shallow)
		{
			$args[] = "-mosaic";
		}
		
		if(db_is_enabled("NGSD"))
		{
			$db = DB::getInstance("NGSD", false);
			$ps_id = get_processed_sample_id($db, $name, false);
			if ($ps_id!=-1)
			{
				$gender = $db->getValue("SELECT s.gender FROM processed_sample ps, sample s WHERE ps.sample_id=s.id AND ps.id={$ps_id}");
				$args[] = "-gender {$gender}";
			}
		}
		
		$parser->execTool("NGS/vc_clincnv_germline.php", implode(" ", $args), true);
		
		//copy results to output folder
		if (file_exists($cnv_out)) $parser->moveFile($cnv_out, $cnvfile);
		if (file_exists($cnv_out2)) $parser->moveFile($cnv_out2, $cnvfile2);
		$mosaic = $folder."/".$name."_mosaic_cnvs.tsv";
		$sample_cnv_name = substr($cnv_out,0,-4);
		$mosaic_out = $sample_cnv_name."_mosaic.tsv";
		if (file_exists($mosaic_out)) $parser->moveFile($mosaic_out, $mosaic);
	}
	else
	{
		check_genome_build($cnvfile, $build);
	}

	// annotate CNV file
	if (file_exists($cnvfile))
	{
		$repository_basedir = repository_basedir();
		$data_folder = get_path("data_folder");
		if ($is_wes || $is_wgs || $is_wgs_shallow) //Genome/Exome: ClinCNV
		{

			$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$cnvfile}", true);
			
		}
		$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$cnvfile}", true);
		$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed -no_duplicates -url_decode -out {$cnvfile}", true);
		$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2024-02.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$cnvfile}", true);


		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2023_3.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$cnvfile}", true);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->exec($ngsbits."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$omim_file} -no_duplicates -url_decode -out {$cnvfile}", true);
		}

		//annotate additional gene info
		$parser->exec($ngsbits."CnvGeneAnnotation", "-in {$cnvfile} -add_simple_gene_names -out {$cnvfile}", true);
		// skip annotation if no connection to the NGSD is possible
		if (db_is_enabled("NGSD"))
		{
			//annotate overlap with pathogenic CNVs
			$parser->exec($ngsbits."NGSDAnnotateCNV", "-in {$cnvfile} -out {$cnvfile}", true);
		}
	}
	else
	{
		trigger_error("CNV file {$cnvfile} does not exist, skipping CNV annotation!", E_USER_WARNING);
	}
}

//structural variants
if (in_array("sv", $steps))
{
	// skip SV calling if only annotation should be done	
	if (!$annotation_only)
	{
		$dragen_output_vcf = $folder."/dragen_variant_calls/{$name}_dragen_svs.vcf.gz";
		if ($use_dragen && get_path("use_dragen_sv_calling"))
		{
			if (!file_exists($dragen_output_vcf)) trigger_error("Dragen SV calling file not found!", E_USER_ERROR);
			if (!in_array("ma", $steps)) trigger_error("'-use_dragen' with no mapping step provided. Using old DRAGEN VCF for SV calling.", E_USER_NOTICE);
					
			//combine BND of INVs to one INV in VCF
			$vcf_inv_corrected = $parser->tempFile("_sv_inv_corrected.vcf");
			$parser->exec(get_path("python27")." ".get_path('manta')."/../libexec/convertInversion.py", get_path("samtools")." {$genome} {$dragen_output_vcf} > {$vcf_inv_corrected}");

			//remove VCF lines with empty "REF". They are sometimes created from convertInversion.py but are not valid
			$vcf_fixed = $parser->tempFile("_sv_fixed.vcf");
			$h = fopen2($vcf_inv_corrected, "r");
			$h2 = fopen2($vcf_fixed, "w");
			while(!feof($h))
			{
				$line = fgets($h);
				$parts = explode("\t", $line);
				if (count($parts)>3 && $parts[3]=="") continue;
				
				fputs($h2, $line);
			}
			fclose($h);
			fclose($h2);

			//sort variants
			$vcf_sorted = $parser->tempFile("_sv_sorted.vcf");
			$parser->exec($ngsbits."VcfSort", "-in {$vcf_fixed} -out {$vcf_sorted}", true);
						
			//mark off-target variants
			if($has_roi)
			{
				$vcf_marked = $parser->tempFile("_sv_marked.vcf");
				$parser->exec($ngsbits."VariantFilterRegions", "-in {$vcf_sorted} -reg ".$sys['target_file']." -mark off-target -out {$vcf_marked}");
				$vcf_sorted = $vcf_marked;
			}
	
			//bgzip and index
			$parser->exec("bgzip", "-c $vcf_marked > $sv_manta_file");
			$parser->exec("tabix", "-p vcf $sv_manta_file", false); //no output logging, because Toolbase::extractVersion() does not return
		}
		else
		{
			//SV calling with manta
			$manta_evidence_dir = "{$folder}/manta_evid";
			create_directory($manta_evidence_dir);

			$manta_args = [
				"-bam ".$used_bam_or_cram,
				"-evid_dir ".$manta_evidence_dir,
				"-out ".$sv_manta_file,
				"-threads ".$threads,
				"-build ".$build,
				"--log ".$parser->getLogFile()
			];
			if($has_roi) $manta_args[] = "-target ".$sys['target_file'];
			if(!$is_wgs) $manta_args[] = "-exome";
			
			$parser->execTool("NGS/vc_manta.php", implode(" ", $manta_args));

			//rename Manta evidence file
			$parser->moveFile("$manta_evidence_dir/evidence_0.$name.bam", "$manta_evidence_dir/{$name}_manta_evidence.bam");
			$parser->moveFile("$manta_evidence_dir/evidence_0.$name.bam.bai", "$manta_evidence_dir/{$name}_manta_evidence.bam.bai");
		}
		
		//create BEDPE files
		$parser->exec("{$ngsbits}VcfToBedpe", "-in $sv_manta_file -out $bedpe_out", true);
	}
	else
	{
		check_genome_build($sv_manta_file, $build);
	}


	//add gene info annotation
	if (db_is_enabled("NGSD"))
	{
		$parser->exec("{$ngsbits}BedpeGeneAnnotation", "-in $bedpe_out -out $bedpe_out -add_simple_gene_names", true);
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
		$parser->exec("{$ngsbits}BedpeAnnotateCounts", "-in $bedpe_out -out $bedpe_out -processing_system ".$sys["name_short"]." -ann_folder {$ngsd_annotation_folder}", true);
		$sys_specific_density_file = $ngsd_annotation_folder."sv_breakpoint_density_".$sys["name_short"].".igv";
		if (file_exists($sys_specific_density_file))
		{
			$parser->exec("{$ngsbits}BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv -density_sys {$sys_specific_density_file}", true);
		}
		else
		{
			$parser->exec("{$ngsbits}BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv", true);
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
		$parser->exec("{$ngsbits}BedpeAnnotateFromBed", "-in $bedpe_out -out $bedpe_out -bed $omim_file -url_decode -replace_underscore -col_name OMIM", true);
	}

	//add CNV overlap annotation
	if (file_exists($cnvfile))
	{
		$parser->exec("{$ngsbits}BedpeAnnotateCnvOverlap", "-in $bedpe_out -out $bedpe_out -cnv $cnvfile", true);
	}

	//update sample entry 
	if (db_is_enabled("NGSD"))
	{
		$bedpe_content = Matrix::fromTSV($bedpe_out);
		$old_comments = $bedpe_content->getComments();
		$new_comments = array();
		foreach ($old_comments as $line) 
		{
			if(starts_with($line, "#SAMPLE="))
			{
				//replace SAMPLE line with updated entry from NGSD
				$new_comments[] = gsvar_sample_header($name, array("DiseaseStatus"=>"Affected"), "#", "\n"); 
			}
			else
			{
				//keep old line
				$new_comments[] = $line;
			}
		}
		$bedpe_content->setComments($new_comments);
		$bedpe_content->toTSV(($bedpe_out));
	}
}

//repeat expansions
if (in_array("re", $steps))
{
	//perform repeat expansion analysis (only for WGS/WES):
	$parser->execTool("NGS/vc_expansionhunter.php", "-in $used_bam_or_cram -out $expansion_hunter_file -build ".$build." -pid $name -threads {$threads}");
}

// create Circos plot - if small variant, CNV or SV calling was done
if ((in_array("vc", $steps) || in_array("cn", $steps) || in_array("sv", $steps)) && !$annotation_only)
{
	if ($is_wes || $is_wgs || $is_wgs_shallow)
	{
		if (file_exists($cnvfile))
		{
			if (file_exists($cnvfile2))
			{
				$parser->execTool("NGS/create_circos_plot.php", "-folder $folder -name $name -build ".$build);
			}
			else
			{
				trigger_error("CNV file $cnvfile2 missing. Cannot create Circos plot!", E_USER_WARNING);
			}
		}
		else
		{
			trigger_error("CNV file $cnvfile missing. Cannot create Circos plot!", E_USER_WARNING);
		}
	}
}

// collect other QC terms
if (in_array("cn", $steps) || in_array("sv", $steps) || in_array("db", $steps))
{
	$terms = [];
	$sources = [];
	
	//CNVs
	if (file_exists($cnvfile))
	{
		$cnv_count_hq = 0;
		$cnv_count_hq_autosomes = 0;
		$cnv_count_loss = 0;
		$cnv_count_gain = 0;
		$h = fopen2($cnvfile, 'r');
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
		
		$sources[] = $cnvfile;
	}

	//SVs
	if (file_exists($bedpe_out))
	{
		$sv_count_pass = 0;
		$sv_count_del = 0;
		$sv_count_dup = 0;
		$sv_count_ins = 0;
		$sv_count_inv = 0;
		$sv_count_bnd = 0;
		$h = fopen2($bedpe_out, 'r');
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
		
		$sources[] = $bedpe_out;
	}
	
	//create qcML file
	if (count($sources)>0)
	{
		$tmp = $parser->tempFile("qc.tsv");
		file_put_contents($tmp, implode("\n", $terms));
		$parser->exec("{$ngsbits}TsvToQC", "-in $tmp -out $qc_other -sources ".implode(" ", $sources));
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
	if(!$somatic) $parser->execTool("NGS/db_check_gender.php", "-in $used_bam_or_cram -pid $name");	
	
	//import variants
	$args = ["-ps {$name}"];
	$import = false;
	if (file_exists($varfile) && !$is_wgs_shallow)
	{
		//check genome build
		check_genome_build($varfile, $build);
		
		$args[] = "-var {$varfile}";
		$args[] = "-var_force";
		$import = true;
	}
	if (file_exists($cnvfile))
	{
		//check genome build
		//this is not possible for CNVs because the file does not contain any information about it
		
		$args[] = "-cnv {$cnvfile}";
		$args[] = "-cnv_force";
		$import = true;
	}
	if (file_exists($bedpe_out))
	{
		//check genome build
		check_genome_build($bedpe_out, $build);
		
		$args[] = "-sv {$bedpe_out}";
		$args[] = "-sv_force";
		$import = true;
	}
	if (file_exists($expansion_hunter_file))
	{
		$args[] = "-re {$expansion_hunter_file}";
		$args[] = "-re_force";
		$import = true;
	}
	if ($import)
	{
		$parser->exec("{$ngsbits}NGSDAddVariantsGermline", implode(" ", $args), true);
	}
}

?>
