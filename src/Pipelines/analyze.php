<?php

/**
	@page analyze
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("analyze", "Complete NGS analysis pipeline.");
$parser->addInfile("folder", "Analysis data folder.", false);
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
$parser->addFlag("no_gender_check", "Skip gender check (done between mapping and variant calling).");
$parser->addFlag("correction_n", "Use Ns for errors by barcode correction.");
$parser->addFlag("somatic", "Set somatic single sample analysis options (i.e. correction_n, clip_overlap).");
$parser->addFlag("annotation_only", "Performs only a reannotation of the already created variant calls.");
$parser->addFlag("no_dragen", "Do not use Illumina DRAGEN calls from NovaSeq X or Dragen server.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
$parser->addFlag("no_splice", "Skip SpliceAI scoring of variants that are not precalculated.");
$parser->addFlag("no_circos", "Skip Circos plot generation.");
$parser->addFlag("no_qc", "Skip calculation of QC metrics.");
$parser->addFlag("no_lowcov", "Skip generation of low-coverage report");
$parser->addString("rna_sample", "Processed sample name of the RNA sample which should be used for annotation.", true, "");
extract($parser->parse($argv));

// create logfile in output folder if no filepath is provided:
if (!file_exists($folder))
{
	exec2("mkdir -p $folder");
}

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($folder."/analyze_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//determine processing system
$sys = load_system($system, $name);
$is_wes = $sys['type']=="WES";
$is_wgs = $sys['type']=="WGS";
$is_panel = $sys['type']=="Panel" || $sys['type']=="Panel Haloplex";
$is_wgs_shallow = $sys['type']=="WGS (shallow)";
$roi = trim($sys['target_file']);
$build = $sys['build'];

//check that ROI is sorted
if ($roi!="")
{
	$roi = realpath($roi);
	if (!bed_is_sorted($roi)) trigger_error("Target region file not sorted: ".$roi, E_USER_ERROR);
}

//disable abra and soft-clipping if DeepVariant is used for calling
$use_freebayes = get_path("use_freebayes");
if (!$use_freebayes)
{
	$no_abra = true;
	$clip_overlap = false;
}

//handle somatic flag
if ($somatic)
{
	$no_abra = true;
	$clip_overlap = true;
	$correction_n = true;
}

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//remove invalid steps
if (in_array("cn", $steps) && $roi=="")
{
	trigger_error("Skipping step 'cn' - Copy number analysis is only supported for processing systems with target region BED file!", E_USER_NOTICE);
	if (($key = array_search("cn", $steps)) !== false) unset($steps[$key]);
}
if (in_array("sv", $steps) && $is_wgs_shallow)
{
	trigger_error("Skipping step 'sv' - Structural variant calling is not supported for shallow WGS samples!", E_USER_NOTICE);
	if (($key = array_search("sv", $steps)) !== false) unset($steps[$key]);
}
if (in_array("re", $steps) && $is_wgs_shallow)
{
	trigger_error("Skipping step 're' - Repeat expansion calling is not supported for shallow WGS samples!", E_USER_NOTICE);
	if (($key = array_search("re", $steps)) !== false) unset($steps[$key]);
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

$genome = genome_fasta($build);

//dragen files
$dragen_folder = "{$folder}/dragen/";
$dragen_bam = "{$dragen_folder}/{$name}.bam";
$dragen_cram = "{$dragen_folder}/{$name}.cram";
$dragen_bam_or_cram_exists = file_exists($dragen_bam) || file_exists($dragen_cram);
$dragen_output_vcf = "{$dragen_folder}/{$name}.hard-filtered.vcf.gz";
$dragen_output_sv_vcf = "{$dragen_folder}/{$name}.sv.vcf.gz";

//output file names:
//mapping
$bamfile = $folder."/".$name.".bam";
$cramfile = $folder."/".$name.".cram";
$bam_or_cram_exists = file_exists($bamfile) || file_exists($cramfile);
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
$sv_manta_file = $folder ."/". $name . "_var_structural_variants.vcf.gz";
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

//prevent accidentally re-mapping if DRAGEN already ran
if (in_array("ma", $steps) && !$no_dragen && file_exists($dragen_folder) && ($bam_or_cram_exists || $dragen_bam_or_cram_exists)) 
{
	trigger_error("'ma' step requested, but sample is already analyzed with DRAGEN. Use '-no_dragen' if you really want to do a re-mapping!", E_USER_ERROR);
}

//move BAM/CRAM from DRAGEN folder to sample folder (on first analysis)
if (!in_array("ma", $steps) && !$no_dragen && file_exists($dragen_folder))
{
	if ($dragen_bam_or_cram_exists)
	{
		if (file_exists($dragen_cram))
		{
			$parser->moveFile($dragen_cram, $cramfile);
			$parser->moveFile($dragen_cram.".crai", $cramfile.".crai");
		}
		else
		{
			$parser->moveFile($dragen_bam, $bamfile);
			$parser->moveFile($dragen_bam.".bai", $bamfile.".bai");
		}
		$bam_or_cram_exists = true;
	}
	else if (!$bam_or_cram_exists)
	{
		trigger_error("Anaylsis without mapping requested, but no BAM/CRAM file found in {$dragen_folder}", E_USER_ERROR);
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
		$parser->execApptainer("ngs-bits", "BamToFastq", "-in $source_file -ref $cram_genome -out1 $tmp1 -out2 $tmp2", [$genome, $folder]);
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
	if($somatic) $args[] = "-somatic_custom_map";
	$used_bam_or_cram = $parser->tempFolder("local_bam")."/".$name.".bam"; //local copy of BAM file to reduce IO over network when mapping is done 
	$parser->execTool("Tools/mapping.php", "-in_for ".implode(" ", $files1)." -in_rev ".implode(" ", $files2)." -system $system -out_folder $folder -out_name $name -local_bam $used_bam_or_cram ".implode(" ", $args)." -threads $threads");
	
	//low-coverage report
	if ($roi!="" && !$is_wgs_shallow && !$no_lowcov)
	{	
		$parser->execApptainer("ngs-bits", "BedLowCoverage", "-in {$roi} -bam $used_bam_or_cram -out $lowcov_file -cutoff 20 -threads {$threads} -ref {$genome}", [$folder, $roi, $genome]);
		if (db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", [$folder]);
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
else if ($bam_or_cram_exists)
{	
	//set BAM/CRAM to use
	$used_bam_or_cram = file_exists($bamfile) ? $bamfile : $cramfile;
	
	//check genome build of BAM
	check_genome_build($used_bam_or_cram, $build);

	//QC for samples mapped/called on NovaSeq X
	if(!file_exists($qc_map) && !$no_qc)
	{
		//QC
		$in_files = array();
		$params = array("-in $used_bam_or_cram", "-out {$qc_map}", "-ref {$genome}", "-build ".ngsbits_build($sys['build']));
		$in_files[] = $folder;
		$in_files[] = $genome;
		if ($roi=="" || $sys['type']=="WGS" || $sys['type']=="WGS (shallow)")
		{
			$params[] = "-wgs";
		}
		else
		{
			$params[] = "-roi {$roi}";
			$in_files[] = $roi;
		}
		if ($sys['build']!="GRCh38")
		{
			$params[] = "-no_cont";
		}
		if ($somatic && file_exists($somatic_custom_panel))
		{
			$params[] = "-somatic_custom_bed $somatic_custom_panel";
			$in_files[] = $somatic_custom_panel;
		}
		if (!file_exists($qc_fastq))
		{
			$params[] = "-read_qc $qc_fastq";
		}
		$parser->execApptainer("ngs-bits", "MappingQC", implode(" ", $params), $in_files);
	}

	//low-coverage regions for samples mapped/called on NovaSeq X / DRAGEN server
	if(!file_exists($lowcov_file) && !$no_lowcov)
	{
		if ($roi!="" && !$is_wgs_shallow)
		{	
			$parser->execApptainer("ngs-bits", "BedLowCoverage", "-in {$roi} -bam $used_bam_or_cram -out $lowcov_file -cutoff 20 -threads {$threads} -ref {$genome}", [$folder, $roi, $genome]);
			if (db_is_enabled("NGSD"))
			{

				$parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", [$folder]);
			}
		}
	}

	//delete fastq files after mapping - remove files only if sample doesn't contain UMIs and the corresponding setting is set in the settings.ini
	if ($sys['umi_type']=="n/a" && get_path("delete_fastq_files")) 
	{
		//check if project overwrites the settings
		$preserve_fastqs = true;
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
			//check for remaining FastQ/ORAs and delete them
			$fastq_files = glob($folder."/*_R?_00?.fastq.gz");
			$ora_files = glob($folder."/*_R?_00?.fastq.ora");

			if ((count($ora_files) + count($fastq_files)) > 0)
			{
				$bam_read_count = get_read_count($used_bam_or_cram, max(8, $threads), array("-F", "2304"), $sys['build']);

				if (count($ora_files)  > 0) $fastq_read_count = get_ora_read_count($ora_files);
				else $fastq_read_count = get_fastq_read_count($fastq_files);

				//calculate relative difference
				$diff = abs($bam_read_count - $fastq_read_count);
				$rel_diff = $diff /(($bam_read_count + $fastq_read_count)/2);

				if (($bam_read_count == $fastq_read_count) || ($rel_diff < 0.001))
				{
					// remove old BAM(s)
					if ($bam_read_count == $fastq_read_count) trigger_error("Read count of FASTQ/ORA files and BAM/CRAM match. Deleting FASTQ/ORA files...", E_USER_NOTICE);
					else if ($rel_diff < 0.001) trigger_error("Read count of FASTQ/ORA files and BAM/CRAM in allowed tolerance (<0.01%) (".($rel_diff*100)."%). Deleting FASTQ/ORA files...", E_USER_NOTICE);
					
					if (count($ora_files)  > 0) $parser->exec("rm", implode(" ", $ora_files));
					else $parser->exec("rm", implode(" ", $fastq_files));
				}
				else
				{
					trigger_error("Cannot delete FASTQ/ORA file(s) - BAM file read counts doesn't match FASTQ/ORA file(s) (".($rel_diff*100)."%).", E_USER_ERROR);
				}
			}
		}
	}
}
else
{
	trigger_error("No BAM/CRAM file found for analysis!", E_USER_ERROR);
}

//check gender after mapping
if(db_is_enabled("NGSD") && $used_bam_or_cram!="" && !$no_gender_check)
{
	
	$db = DB::getInstance("NGSD", false);
	$info = get_processed_sample_info($db, $name, false);
	
	if (!is_null($info) && $info["is_tumor"] == "0")
	{
		$args = array();
		$args[] = "-in $used_bam_or_cram";
		$args[] = "-pid $name";
		
		//extra check for germline ffpe samples, as they often have high SNV deviation
		if ($info["is_ffpe"] == "1")
		{
			$args[] = "-check_sry_cov";
		}
		
		$parser->execTool("Tools/db_check_gender.php", implode(" ", $args));
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
		if ($roi!="") 
		{
			$only_mito_in_target_region = exec2("cat {$roi} | cut -f1 | uniq")[0][0] == "chrMT";
		}
		//activate mito-calling for sWGS
		if ($is_wgs_shallow) $only_mito_in_target_region = true;

		//perform main variant calling on autosomes/genosomes
		if(!$only_mito_in_target_region)
		{			
			if (!$no_dragen && file_exists($dragen_output_vcf))
			{
				trigger_error("DRAGEN analysis found in sample folder. Using this data for small variant calling. ", E_USER_NOTICE);
				$pipeline = [];

				$pipeline[] = array("zcat", $dragen_output_vcf);
				
				//filter by target region (extended by 200) and quality 5
				$target = $parser->tempFile("_roi_extended.bed");
				$parser->execApptainer("ngs-bits", "BedExtend", "-in {$roi} -n 200 -out $target -fai ".$genome.".fai", [$roi, $genome]);
				$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfFilter", "-reg {$target} -qual 5 -filter_clear -remove_invalid -ref $genome", [$genome], [], true));

				//split multi-allelic variants
				$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "", [], [], true)];

				//normalize all variants and align INDELs to the left
				$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref $genome", [$genome], [], true));

				//sort variants by genomic position
				$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfStreamSort", "", [], [], true));

				//zip
				$pipeline[] = array("", $parser->execApptainer("htslib", "bgzip", "-c > $vcffile", [], [dirname($vcffile)], true));

				//execute pipeline
				$parser->execPipeline($pipeline, "Dragen small variants post processing");

				//mark off-target variants
				$tmp = $parser->tempFile("_offtarget.vcf");
				$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in $vcffile -mark off-target -reg {$roi} -out $tmp", [$folder, $roi]);
				
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
				$parser->execApptainer("htslib", "bgzip", "-c $tmp2 > $vcffile", [], [dirname($vcffile)]);

				//index output file
				$parser->execApptainer("htslib", "tabix", "-p vcf $vcffile", [], [dirname($vcffile)]);
			}
			elseif ($use_freebayes) //perform variant calling with freebayes if set in settings.ini
			{
				$args = [];
				$args[] = "-bam ".$used_bam_or_cram;
				$args[] = "-out ".$vcffile;
				$args[] = "-build ".$build;
				$args[] = "-threads ".$threads;
				if ($roi!="")
				{
					$args[] = "-target {$roi}";
					$args[] = "-target_extend 200";
				}
				$args[] = "-min_af ".$min_af;
				$args[] = "-min_mq ".$min_mq;
				$args[] = "-min_bq ".$min_bq;
				$parser->execTool("Tools/vc_freebayes.php", implode(" ", $args));
			}
			else //perform variant calling with DeepVariant
			{
				$args = [];
				if ($is_wes || $is_panel) $args[] = "-model_type WES";
				else if ($is_wgs) $args[] = "-model_type WGS";
				else trigger_error("Unsupported system type '".$sys['type']."' detected in $system. Compatible system types are: WES, WGS, Panel, Panel Haloplex.", E_USER_ERROR);
				$args[] = "-bam ".$used_bam_or_cram;
				$args[] = "-out ".$vcffile;
				$args[] = "-build ".$build;
				$args[] = "-threads ".$threads;
				if ($roi!="")
				{
					$args[] = "-target {$roi}";
					$args[] = "-target_extend 200";
				}
				$args[] = "-min_af ".$min_af;
				$args[] = "-min_mq ".$min_mq;
				$args[] = "-min_bq ".$min_bq;
				$args[] = "-allow_empty_examples";

				$parser->execTool("Tools/vc_deepvariant.php", implode(" ", $args));
			}
		}
		
		//perform special variant calling for mitochondria
		$mito = enable_special_mito_vc($sys) || $only_mito_in_target_region;
		if ($mito)
		{
			$vcffile_mito = $parser->tempFile("_mito.vcf.gz");
			if (!$no_dragen && file_exists($dragen_output_vcf) && contains_mito($dragen_output_vcf))
			{
				trigger_error("DRAGEN analysis found in sample folder. Using this data for mito small variant calling. ", E_USER_NOTICE);
				$pipeline = [];
				
				$pipeline[] = array("zcat", $dragen_output_vcf);
				
				//filter by target region and quality 5
				$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfFilter", "-reg chrMT:1-16569 -qual 5 -filter_clear -ref $genome", [$genome], [], true));

				//split multi-allelic variants
				$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "", [], [], true)];

				//normalize all variants and align INDELs to the left
				$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref $genome", [$genome], [], true));

				//sort variants by genomic position
				$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfStreamSort", "", [], [], true));

				//zip
				$pipeline[] = array("", $parser->execApptainer("htslib", "bgzip", "-c > $vcffile_mito", [], [], true));

				//execute pipeline
				$parser->execPipeline($pipeline, "Dragen small variants post processing (chrMT)");

				//index output file
				$parser->execApptainer("htslib", "tabix", "-p vcf $vcffile_mito");
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
				$parser->execTool("Tools/vc_mosaic.php", implode(" ", $args));
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
		if ($is_wgs || ($is_wes && $roi!="") || ($is_panel && $roi!=""))
		{
			//determine region
			if ($is_wgs)
			{
				$roi_low_mappabilty = repository_basedir()."data/misc/low_mappability_region/wes_mapq_eq0.bed";
			}
			else
			{
				$roi_low_mappabilty = $parser->tempFile("_mapq0.bed");
				$parser->execApptainer("ngs-bits", "BedIntersect", "-in {$roi} -in2 ".repository_basedir()."data/misc/low_mappability_region/wes_mapq_eq0.bed -out $roi_low_mappabilty", [$roi, repository_basedir()."data/misc/low_mappability_region/wes_mapq_eq0.bed"]);
			}
			
			if (bed_size($roi_low_mappabilty)>0)
			{
				$tmp_low_mappability = $parser->tempFile("_low_mappability.vcf.gz");

				if ($use_freebayes) //perform variant calling with freebayes if set in settings.ini
				{
					$args = [];
					$args[] = "-bam ".$used_bam_or_cram;
					$args[] = "-out ".$tmp_low_mappability;
					$args[] = "-build ".$build;
					$args[] = "-threads ".$threads;
					$args[] = "-target ".$roi_low_mappabilty;
					$args[] = "-min_af ".$min_af;
					$args[] = "-min_mq 0";
					$args[] = "-min_bq ".$min_bq;
					$parser->execTool("Tools/vc_freebayes.php", implode(" ", $args));
				}
				else
				{
					$args = [];
					if ($is_wes || $is_panel) $args[] = "-model_type WES";
					else if ($is_wgs) $args[] = "-model_type WGS";
					else trigger_error("Unsupported system type '".$sys['type']."' detected in {$system}. Compatible system types are: WES, WGS, Panel, Panel Haloplex.", E_USER_ERROR);
					$args[] = "-bam ".$used_bam_or_cram;
					$args[] = "-out ".$tmp_low_mappability;
					$args[] = "-build ".$build;
					$args[] = "-threads ".$threads;
					$args[] = "-target $roi_low_mappabilty";
					$args[] = "-min_af ".$min_af;
					$args[] = "-min_mq 0";
					$args[] = "-min_bq ".$min_bq;
					$args[] = "-allow_empty_examples";

					$parser->execTool("Tools/vc_deepvariant.php", implode(" ", $args));
				}

				//unzip
				$tmp_low_mappability2 = $parser->tempFile("_low_mappability.vcf");
				$parser->exec("zcat", "$tmp_low_mappability > $tmp_low_mappability2", true);
				
				//add to main variant list
				$tmp2 = $parser->tempFile("_merged_low_mappability.vcf");
				$parser->execApptainer("ngs-bits", "VcfAdd", "-in $vcf $tmp_low_mappability2 -skip_duplicates -filter low_mappability -filter_desc Variants_in_reads_with_low_mapping_score. -out $tmp2");
				$parser->moveFile($tmp2, $vcf);
			}
		}
			
		//call mosaic variants on target region (WES or Panel) or exonic/splicing region (WGS)
		if (ngsbits_build($build)!="non_human" && ($is_wgs || ($is_wes && $roi!="") || ($is_panel && $roi!="")))
		{
			$tmp_mosaic = $parser->tempFile("_mosaic.vcf");
			$args = [];
			$args[] = "-in {$used_bam_or_cram}";
			$args[] = "-out {$tmp_mosaic}";
			$args[] = "-no_zip";
			$args[] = "-target ".(($is_panel || $is_wes) ? $roi : repository_basedir()."data/gene_lists/gene_exons_pad20.bed");
			$args[] = "-threads ".$threads;
			$args[] = "-build ".$build;
			$args[] = "-min_af 0.03";
			$args[] = "-min_obs ".($is_wgs ?  "1" : "2");	
			$parser->execTool("Tools/vc_mosaic.php", implode(" ", $args));
			
			//add to main variant list
			$tmp2 = $parser->tempFile("_merged_mosaic.vcf");
			$parser->execApptainer("ngs-bits", "VcfAdd", "-in $vcf $tmp_mosaic -skip_duplicates -filter mosaic -filter_desc Putative_mosaic_variants. -out $tmp2");
			$parser->moveFile($tmp2, $vcf);
		}
				
		//sort and remove unused contig lines
		$parser->execApptainer("ngs-bits", "VcfSort", "-in $vcf -remove_unused_contigs -out $vcf");
		
		//zip and index
		$parser->execApptainer("htslib", "bgzip", "-c $vcf > $vcffile", [], [dirname($vcffile)]);
		$parser->execApptainer("htslib", "tabix", "-f -p vcf $vcffile", [], [dirname($vcffile)]);
		
		//create b-allele frequency file
		$params = array();
		$params[] = "-vcf {$vcffile}";
		$params[] = "-name {$name}";
		$params[] = "-out {$baffile}";
		if ($is_wgs)
		{
			$params[] = "-downsample 100";
		}
		$parser->execTool("Tools/baf_germline.php", implode(" ", $params));
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
	if ($no_splice) $args[] = "-no_splice";
	$parser->execTool("Tools/annotate.php", implode(" ", $args));

	//check for truncated output
	if ($is_wgs) $parser->execTool("Tools/check_for_missing_chromosomes.php", "-in {$vcffile_annotated} -max_missing_perc 5");
		
	//ROH detection
	if ($is_wes || $is_wgs)
	{
		$in_files = [];
		$in_files[] = $folder;
		$in_files[] = repository_basedir()."/data/gene_lists/genes.bed";
		$in_files[] = repository_basedir()."/data/misc/roh_exclude_regions.bed";
		$args = [];
		$args[] = "-in $vcffile_annotated";
		$args[] = "-out $rohfile";
		$args[] = "-var_af_keys gnomADg_AF";
		$args[] = "-exclude ".repository_basedir()."/data/misc/roh_exclude_regions.bed";
		$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file)) $in_files[] = $omim_file;
		$args[] = "-annotate ".repository_basedir()."/data/gene_lists/genes.bed ".(file_exists($omim_file) ? $omim_file : "");
		$parser->execApptainer("ngs-bits", "RohHunter", implode(" ", $args), $in_files);
	}

	//PRS calculation 
	if ($is_wgs)
	{
		$prs_folder = repository_basedir()."/data/misc/prs/";
		$prs_scoring_files = glob($prs_folder."/*_".$build.".vcf");
		if (count($prs_scoring_files) > 0)
		{
			$parser->execApptainer("ngs-bits", "VcfCalculatePRS", "-in $vcffile -bam $used_bam_or_cram -out $prsfile -prs ".implode(" ", $prs_scoring_files)." -ref $genome", [$folder, $genome, $prs_folder]);
		}
	}
	
	//determine ancestry
	if (ngsbits_build($build) != "non_human")
	{
		$parser->execApptainer("ngs-bits", "SampleAncestry", "-in {$vcffile} -out {$ancestry_file} -build ".ngsbits_build($build), [$folder]);
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
						["", $parser->execApptainer("ngs-bits", "BedChunk", "-in {$roi} -n {$bin_size}", [$roi], [], true)],
						["", $parser->execApptainer("ngs-bits", "BedAnnotateGC", "-clear -ref ".$genome, [$genome], [], true)],
						["", $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-out {$bed}", [], [dirname($bed)], true)]
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
						["", $parser->execApptainer("ngs-bits", "BedAnnotateGC", "-in {$roi} -clear -ref ".$genome, [$roi, $genome], [], true)],
						["", $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-out {$bed}", [], [dirname($bed)], true)]
					];
				$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
			}
		}

		//create coverage profile
		$tmp_folder = $parser->tempFolder();
		$cov_file = $cov_folder."/{$name}.cov.gz";
		$cov_tmp_unzipped = $tmp_folder."/{$name}.cov";
		$cov_tmp = $cov_tmp_unzipped.".gz";
		$parser->execApptainer("ngs-bits", "BedCoverage", "-clear -min_mapq 0 -decimals 4 -bam {$used_bam_or_cram} -in {$bed} -out {$cov_tmp_unzipped} -threads {$threads} -ref {$genome}", [$genome, $bed, $folder]);
		$parser->exec("gzip", "-9 {$cov_tmp_unzipped}");
		
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
		
		$parser->execTool("Tools/vc_clincnv_germline.php", implode(" ", $args), true);
		
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

			$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnvfile} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$cnvfile}", [$folder, "{$repository_basedir}/data/misc/af_genomes_imgag.bed"]);
			
		}
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnvfile} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$cnvfile}", [$folder, "{$repository_basedir}/data/misc/cn_pathogenic.bed"]);
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnvfile} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed -no_duplicates -url_decode -out {$cnvfile}", [$folder, "{$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed"]);
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnvfile} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2025-09.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$cnvfile}", [$folder, "{$data_folder}/dbs/ClinVar/clinvar_cnvs_2025-09.bed"]);


		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2025_2.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnvfile} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$cnvfile}", [$folder, $hgmd_file]);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnvfile} -in2 {$omim_file} -no_duplicates -url_decode -out {$cnvfile}", [$folder, $omim_file]);
		}

		//annotate gene info
		if(db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "CnvGeneAnnotation", "-in {$cnvfile} -add_simple_gene_names -out {$cnvfile}", [$folder]);
		}
		
		//annotate overlap with pathogenic CNVs
		if (db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "NGSDAnnotateCNV", "-in {$cnvfile} -out {$cnvfile}", [$folder]);
		}
		
		//check for truncated output
		if ($is_wgs) $parser->execTool("Tools/check_for_missing_chromosomes.php", "-in {$cnvfile}");
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
		
		if (!$no_dragen && file_exists($dragen_output_sv_vcf))
		{
			trigger_error("DRAGEN analysis found in sample folder. Using this data for structural variant calling. ", E_USER_NOTICE);
					
			//combine BND of INVs to one INV in VCF
			$vcf_inv_corrected = $parser->tempFile("_sv_inv_corrected.vcf");
			$inv_script = repository_basedir()."/src/Tools/convertInversion.py";
			$vc_manta_command = "python2 ".$inv_script;
			$vc_manta_parameters = "/usr/local/bin/samtools {$genome} {$dragen_output_sv_vcf} dragen > {$vcf_inv_corrected}";
			$parser->execApptainer("manta", $vc_manta_command, $vc_manta_parameters, [$genome, $inv_script, $dragen_output_sv_vcf]);

			// fix VCF file (remove variants with empty "REF" entry and duplicates)
			$vcf_fixed = $parser->tempFile("_sv_fixed.vcf");
			$parser->execApptainer("ngs-bits", "MantaVcfFix", "-in {$vcf_inv_corrected} -out {$vcf_fixed}");

			//sort variants
			$vcf_sorted = $parser->tempFile("_sv_sorted.vcf");
			$parser->execApptainer("ngs-bits", "VcfSort", "-in {$vcf_fixed} -out {$vcf_sorted}");
						
			//mark off-target variants
			if($roi!="")
			{
				$vcf_marked = $parser->tempFile("_sv_marked.vcf");
				$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in {$vcf_sorted} -reg {$roi} -mark off-target -out {$vcf_marked}", [$roi]);
				$vcf_sorted = $vcf_marked;
			}
	
			//bgzip and index
			$parser->execApptainer("htslib", "bgzip", "-c $vcf_marked > $sv_manta_file", [], [dirname($sv_manta_file)]);
			$parser->execApptainer("htslib", "tabix", "-p vcf $sv_manta_file", [], [dirname($sv_manta_file)]);
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
			if($roi!="") $manta_args[] = "-target {$roi}";
			if(!$is_wgs) $manta_args[] = "-exome";
			
			$parser->execTool("Tools/vc_manta.php", implode(" ", $manta_args));

			//rename Manta evidence file
			$parser->moveFile("$manta_evidence_dir/evidence_0.$name.bam", "$manta_evidence_dir/{$name}_manta_evidence.bam");
			$parser->moveFile("$manta_evidence_dir/evidence_0.$name.bam.bai", "$manta_evidence_dir/{$name}_manta_evidence.bam.bai");
		}
		
		//create BEDPE files
		$parser->execApptainer("ngs-bits", "VcfToBedpe", "-in $sv_manta_file -out $bedpe_out", [$folder]);
	}
	else
	{
		check_genome_build($sv_manta_file, $build);
	}


	//add gene info annotation
	if (db_is_enabled("NGSD"))
	{
		$parser->execApptainer("ngs-bits", "BedpeGeneAnnotation", "-in $bedpe_out -out $bedpe_out -add_simple_gene_names", [$folder]);
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
			"-in $bedpe_out",
			"-out $bedpe_out",
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
			$parser->execApptainer("ngs-bits", "BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv -density_sys {$sys_specific_density_file}", [$folder, $ngsd_annotation_folder]);
		}
		else
		{
			$parser->execApptainer("ngs-bits", "BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv", [$folder, $ngsd_annotation_folder]);
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
			$parser->execApptainer("ngs-bits", "NGSDAnnotateSV", "-in {$bedpe_out} -out {$bedpe_out}", [$folder]);
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
		$parser->execApptainer("ngs-bits", "BedpeAnnotateFromBed", "-in $bedpe_out -out $bedpe_out -bed $omim_file -url_decode -replace_underscore -col_name OMIM", [$folder, $omim_file]);
	}

	//add CNV overlap annotation
	if (file_exists($cnvfile))
	{
		$parser->execApptainer("ngs-bits", "BedpeAnnotateCnvOverlap", "-in $bedpe_out -out $bedpe_out -cnv $cnvfile", [$folder]);
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
	
	//check for truncated output
	if ($is_wgs) $parser->execTool("Tools/check_for_missing_chromosomes.php", "-in {$bedpe_out}");
}

//repeat expansions
if (in_array("re", $steps))
{
	//perform repeat expansion analysis (only for WGS/WES):
	$parser->execTool("Tools/vc_expansionhunter.php", "-in $used_bam_or_cram -out $expansion_hunter_file -build ".$build." -pid $name -threads {$threads}");
}

// create Circos plot - if small variant, CNV or SV calling was done
if ((in_array("vc", $steps) || in_array("cn", $steps) || in_array("sv", $steps)) && !$annotation_only && !$no_circos)
{
	if ($is_wes || $is_wgs || $is_wgs_shallow)
	{
		if (file_exists($cnvfile))
		{
			if (file_exists($cnvfile2))
			{
				$parser->execTool("Tools/create_circos_plot.php", "-folder $folder -name $name -build ".$build);
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
if ((in_array("cn", $steps) || in_array("sv", $steps) || in_array("db", $steps)) && !$no_qc)
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
		$parser->execApptainer("ngs-bits", "TsvToQC", "-in $tmp -out $qc_other -sources ".implode(" ", $sources), [$folder]);
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
	
	//import variants
	$args = [];
	if (file_exists($varfile) && !$is_wgs_shallow)
	{
		//check genome build
		check_genome_build($varfile, $build);
		
		$args[] = "-var {$varfile}";
	}
	if (file_exists($cnvfile))
	{
		//check genome build
		//this is not possible for CNVs because the file does not contain any information about it
		
		$args[] = "-cnv {$cnvfile}";
	}
	if (file_exists($bedpe_out))
	{
		//check genome build
		check_genome_build($bedpe_out, $build);
		
		$args[] = "-sv {$bedpe_out}";
	}
	if (file_exists($expansion_hunter_file))
	{
		$args[] = "-re {$expansion_hunter_file}";
	}
	if (count($args)>0)
	{
		$args[] = "-ps {$name}";
		$parser->execApptainer("ngs-bits", "NGSDAddVariantsGermline", implode(" ", $args), [$folder]);
	}
}

?>
