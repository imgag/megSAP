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
$steps_all = array("ma", "vc", "cn", "sv", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, cn=copy-number analysis, sv=structural-variant analysis, db=import into NGSD.", true, "ma,vc,cn,sv,db");
$parser->addFloat("min_af", "Minimum VAF cutoff used for variant calling (freebayes 'min-alternate-fraction' parameter).", true, 0.1);
$parser->addFloat("min_bq", "Minimum base quality used for variant calling (freebayes 'min-base-quality' parameter).", true, 15);
$parser->addFloat("min_mq", "Minimum mapping quality used for variant calling (freebayes 'min-mapping-quality' parameter).", true, 1);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("clip_overlap", "Soft-clip overlapping read pairs.");
$parser->addFlag("no_abra", "Skip realignment with ABRA.");
$parser->addFlag("no_trim", "Skip adapter trimming with SeqPurge.");
$parser->addFlag("start_with_abra", "Skip all steps before indel realignment of BAM file.");
$parser->addFlag("correction_n", "Use Ns for errors by barcode correction.");
$parser->addFlag("somatic", "Set somatic single sample analysis options (i.e. correction_n, clip_overlap).");
$parser->addFlag("annotation_only", "Performs only a reannotation of the already created variant calls.");
$parser->addFlag("use_dragen", "Use Illumina DRAGEN server for mapping instead of standard BWA-MEM.");
extract($parser->parse($argv));

// create logfile in output folder if no filepath is provided:
if ($parser->getLogFile() == "") $parser->setLogFile($folder."/analyze_".date("YmdHis").".log");

//init
if($folder=="default")
{
	$folder = $folder;
}
$ngsbits = get_path("ngs-bits");

//determine processing system
$sys = load_system($system, $name);
$is_wes = $sys['type']=="WES";
$is_wgs = $sys['type']=="WGS";
$is_wgs_shallow = $sys['type']=="WGS (shallow)";
$has_roi = $sys['target_file']!="";

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

//check user for DRAGEN mapping
if ($use_dragen && ($user != get_path("dragen_user")))
{
	trigger_error("Analysis has to be run as user '".get_path("dragen_user")."' if DRAGEN mapping should be used!", E_USER_ERROR);
}

//remove invalid steps
if (in_array("vc", $steps) && $is_wgs_shallow)
{
	trigger_error("Skipping step 'vc' - Variant calling is not supported for shallow WGS samples!", E_USER_NOTICE);
	if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
}

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
$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);


//output file names:
//mapping
$bamfile = $folder."/".$name.".bam";
$lowcov_file = $folder."/".$name."_".$sys["name_short"]."_lowcov.bed";
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


// folder for local BAM file
$bam_folder = $parser->tempFolder("local_bam");
$local_bamfile = $bam_folder."/".$name.".bam";

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


// copy BAM file to tmp for all calling steps
if (!$annotation_only && !in_array("ma", $steps))
{
	if(in_array("vc", $steps) || in_array("cn", $steps) || in_array("sv", $steps))
	{
		if (!file_exists($bamfile)) trigger_error("No BAM file found in Sample folder! Cannot perform any calling steps!", E_USER_ERROR);
		if (!file_exists($bamfile.".bai")) trigger_error("No BAM index file found in Sample folder!", E_USER_ERROR);
		// copy BAM file to local tmp
		$parser->copyFile($bamfile, $local_bamfile);
		$parser->copyFile($bamfile.".bai", $local_bamfile.".bai");
	}
}

//mapping
if (in_array("ma", $steps))
{
	//determine input FASTQ files
	$in_for = $folder."/*_R1_00?.fastq.gz";
	$in_rev = $folder."/*_R2_00?.fastq.gz";
	$in_index = $folder."/*_index_*.fastq.gz";
	
	//find FastQ input files
	$files1 = glob($in_for);
	$files2 = glob($in_rev);
	$files_index = glob($in_index);
	if (count($files1)!=count($files2))
	{
		trigger_error("Found mismatching forward and reverse read file count!\n Forward: $in_for\n Reverse: $in_rev.", E_USER_ERROR);
	}
	if (!$start_with_abra)
	{
		if (count($files1)==0)
		{
			if(file_exists($bamfile))
			{
				trigger_error("No FASTQ files found in folder. Using BAM file to generate FASTQ files.", E_USER_NOTICE);

				// extract reads from BAM file
				$in_fq_for = $folder."/{$name}_BamToFastq_R1_001.fastq.gz";
				$in_fq_rev = $folder."/{$name}_BamToFastq_R2_001.fastq.gz";
				$tmp1 = $parser->tempFile(".fastq.gz");
				$tmp2 = $parser->tempFile(".fastq.gz");
				$parser->exec("{$ngsbits}BamToFastq", "-in $bamfile -out1 $tmp1 -out2 $tmp2", true);
				$parser->moveFile($tmp1, $in_fq_for);
				$parser->moveFile($tmp2, $in_fq_rev);
				
				// use generated fastq files for mapping
				$files1 = array($in_fq_for);
				$files2 = array($in_fq_rev);
			}
			else
			{
				trigger_error("Found no read files found matching '$in_for' or '$in_rev'!", E_USER_ERROR);
			}
		}
	}
	
	$args = array();
	if($clip_overlap) $args[] = "-clip_overlap";
	if($no_abra) $args[] = "-no_abra";
	if($no_trim) $args[] = "-no_trim";
	if($start_with_abra) $args[] = "-start_with_abra";
	if($correction_n) $args[] = "-correction_n";
	if(!empty($files_index)) $args[] = "-in_index " . implode(" ", $files_index);
	if($use_dragen) $args[] = "-use_dragen";
	$parser->execTool("Pipelines/mapping.php", "-in_for ".implode(" ", $files1)." -in_rev ".implode(" ", $files2)." -system $system -out_folder $folder -out_name $name -local_bam $local_bamfile ".implode(" ", $args)." -threads $threads");

	//low-coverage report
	if ($has_roi && !$is_wgs && !$is_wgs_shallow)
	{	
		$parser->exec("{$ngsbits}BedLowCoverage", "-in ".$sys['target_file']." -bam $local_bamfile -out $lowcov_file -cutoff 20", true);
		if (db_is_enabled("NGSD"))
		{
			$parser->exec("{$ngsbits}BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", true);
		}
	}

	//delete fastq files after mapping
	$delete_fastq_files = get_path("delete_fastq_files", true);
	// remove files only if sample doesn't contain UMIs and the corresponding setting is set in the settings.ini
	if (($sys['umi_type'] == "n/a") && ($delete_fastq_files==true || $delete_fastq_files=="true")) 
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
			//check if BAM and BAM index exists:
			$bam_exists = file_exists($bamfile) && file_exists($bamfile.".bai"); 
			if ($bam_exists && count($files1)>0 && count($files1)>0)
			{
				//check file sizes:
				//FASTQ
				$fastq_file_size = 0;
				foreach($fastq_files as $fq_file)
				{
					$fastq_file_size += filesize($fq_file);
				}

				//BAM
				$bamfile_size = filesize($bamfile);

				if ($bamfile_size / $fastq_file_size > 0.4)
				{
					// BAM exists and has a propper size: FASTQ files can be deleted
					foreach($fastq_files as $fq_file)
					{
						unlink($fq_file);
					}
				}
				else
				{
					trigger_error("BAM file smaller than 40% of FASTQ", E_USER_ERROR);
				}
			}
			else
			{
				trigger_error("No BAM/BAI file found!", E_USER_ERROR);
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
		$args = array();
		if ($has_roi)
		{
			$args[] = "-target ".$sys['target_file'];
			$args[] = "-target_extend 50";
		}
		$args[] = "-min_af ".$min_af;
		$args[] = "-min_mq ".$min_mq;
		$args[] = "-min_bq ".$min_bq;
		
		//Do not call standard pipeline if there is only mitochondiral chrMT in target region
		$only_mito_in_target_region = false;
		if ($has_roi) $only_mito_in_target_region = exec2("cat ".$sys['target_file']." | cut -f1 | uniq")[0][0] == "chrMT";
		if(!$only_mito_in_target_region)
		{
			$parser->execTool("NGS/vc_freebayes.php", "-bam $local_bamfile -out $vcffile -build ".$sys['build']." -threads $threads ".implode(" ", $args));
		}
		
		//perform special variant calling for mitochondria
		$mito = enable_special_mito_vc($sys) || $only_mito_in_target_region;
		if ($mito)
		{
			$target_mito = $parser->tempFile("_mito.bed");
			file_put_contents($target_mito, "chrMT\t0\t16569");
			
			$args = array();
			$args[] = "-no_ploidy";
			$args[] = "-min_af 0.01";
			$args[] = "-min_mq ".$min_mq;
			$args[] = "-min_bq ".$min_bq;
			$args[] = "-target $target_mito";
			$vcffile_mito = $parser->tempFile("_mito.vcf.gz");
			$parser->execTool("NGS/vc_freebayes.php", "-bam $local_bamfile -out $vcffile_mito -build ".$sys['build']." ".implode(" ", $args));
		}
		
		if($only_mito_in_target_region) 
		{
			$parser->copyFile($vcffile_mito, $vcffile);
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
		
		$parser->exec("bgzip", "-c $vcf > $vcffile", false); //no output logging, because Toolbase::extractVersion() does not return
		$parser->exec("tabix", "-f -p vcf $vcffile", false); //no output logging, because Toolbase::extractVersion() does not return
		
		
		//create b-allele frequency file
		$params = array();
		$params[] = "-vcf {$vcffile}";
		$params[] = "-bam {$local_bamfile}";
		$params[] = "-out {$baffile}";
		$params[] = "-build ".$sys['build'];
		if ($is_wgs)
		{
			$params[] = "-downsample 100";
		}
		$parser->execTool("NGS/baf_germline.php", implode(" ", $params));
	}

	// annotation
	$args = array("-out_name $name", "-out_folder {$folder}", "-system {$system}", "--log ".$parser->getLogFile());
	if ($is_wgs) $args[] = "-updown";
	$args[] = "-threads {$threads}";
	$parser->execTool("Pipelines/annotate.php", implode(" ", $args));
	
	//ROH detection
	if ($is_wes || $is_wgs)
	{
		$args = array();
		$args[] = "-in $vcffile_annotated";
		$args[] = "-out $rohfile";
		$args[] = "-var_af_keys_vep AF,gnomAD_AF -var_af_keys gnomADg_AF"; //use 1000g, gnomAD exome, gnomAD genome
		$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; //optional because of license
		$args[] = "-annotate ".repository_basedir()."/data/gene_lists/genes.bed ".(file_exists($omim_file) ? $omim_file : "");
		$parser->exec("{$ngsbits}RohHunter", implode(" ", $args), true);
	}

	//PRS calculation 
	if ($is_wgs)
	{
		$prs_folder = repository_basedir()."/data/misc/prs/";
		$prs_scoring_files = glob($prs_folder."/*.vcf");
		$parser->exec("{$ngsbits}VcfCalculatePRS", "-in $vcffile -out $prsfile -prs ".implode(" ", $prs_scoring_files), true);
	}
	
	//determine ancestry
	$parser->exec(get_path("ngs-bits")."SampleAncestry", "-in {$vcffile} -out {$ancestry_file}", true);
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
						["{$ngsbits}BedAnnotateGC", "-ref ".genome_fasta($sys['build'])],
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
						["{$ngsbits}BedAnnotateGC", "-in ".$sys['target_file']." -ref ".genome_fasta($sys['build'])],
						["{$ngsbits}BedAnnotateGenes", "-out {$bed}"],
					];
				$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
			}
		}

		//create coverage profile
		$tmp_folder = $parser->tempFolder();
		$cov_file = $cov_folder."/{$name}.cov";
		if (!file_exists($cov_file))
		{

			$parser->log("Calculating coverage file for CN calling...");
			$cov_tmp = $tmp_folder."/{$name}.cov";
			$parser->exec("{$ngsbits}BedCoverage", "-min_mapq 0 -decimals 4 -bam {$local_bamfile} -in {$bed} -out {$cov_tmp}", true);
			
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

		$parser->execTool("NGS/vc_clincnv_germline.php", implode(" ", $args), true);
		
		//copy results to output folder
		if (file_exists($cnv_out)) $parser->moveFile($cnv_out, $cnvfile);
		if (file_exists($cnv_out2)) $parser->moveFile($cnv_out2, $cnvfile2);
		$mosaic = $folder."/".$name."_mosaic_cnvs.tsv";
		$sample_cnv_name = substr($cnv_out,0,-4);
		$mosaic_out = $sample_cnv_name."_mosaic.tsv";
		if (file_exists($mosaic_out)) $parser->moveFile($mosaic_out, $mosaic);

		//create dummy GSvar file for shallow WGS (needed to be able to open the sample in GSvar)
		if ($is_wgs_shallow)
		{
			$content = array(
				"##ANALYSISTYPE=GERMLINE_SINGLESAMPLE",
				"##SAMPLE=<ID={$name},Gender=n/a,ExternalSampleName=n/a,IsTumor=n/a,IsFFPE=n/a,DiseaseGroup=n/a,DiseaseStatus=affected>",
				"##DESCRIPTION={$name}=Genotype of variant in sample.",
				"##DESCRIPTION=filter=Annotations for filtering and ranking variants.",
				"##DESCRIPTION=quality=Quality parameters - Quality parameters - variant quality (QUAL), depth (DP), allele frequency (AF), mean mapping quality of alternate allele (MQM), probability of strand bias for alternate bases as phred score (SAP), probability of allele ballance as phred score (ABP).",
				"##DESCRIPTION=gene=Affected gene list (comma-separated).",
				"##DESCRIPTION=variant_type=Variant type.",
				"##DESCRIPTION=coding_and_splicing=Coding and splicing details (Gene, ENST number, type, impact, exon/intron number, HGVS.c, HGVS.p, Pfam domain).",
				"##DESCRIPTION=regulatory=Regulatory consequence details.",
				"##DESCRIPTION=OMIM=OMIM database annotation.",
				"##DESCRIPTION=ClinVar=ClinVar database annotation.",
				"##DESCRIPTION=HGMD=HGMD database annotation.",
				"##DESCRIPTION=RepeatMasker=RepeatMasker annotation.",
				"##DESCRIPTION=dbSNP=Identifier in dbSNP database.",
				"##DESCRIPTION=1000g=Allele frequency in 1000 genomes project.",
				"##DESCRIPTION=gnomAD=Allele frequency in gnomAD project.",
				"##DESCRIPTION=gnomAD_hom_hemi=Homoyzgous counts and hemizygous counts of gnomAD project (genome data).",
				"##DESCRIPTION=gnomAD_sub=Sub-population allele frequenciens (AFR,AMR,EAS,NFE,SAS) in gnomAD project.",
				"##DESCRIPTION=phyloP=phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6.",
				"##DESCRIPTION=Sift=Sift effect prediction and score for each transcript: D=damaging, T=tolerated.",
				"##DESCRIPTION=PolyPhen=PolyPhen (humVar) effect prediction and score for each transcript: D=probably damaging, P=possibly damaging, B=benign.",
				"##DESCRIPTION=fathmm-MKL=fathmm-MKL score (for coding/non-coding regions). Deleterious threshold > 0.5.",
				"##DESCRIPTION=CADD=CADD pathogenicity prediction scores (scaled phred-like). Deleterious threshold > 10-20.",
				"##DESCRIPTION=REVEL=REVEL pathogenicity prediction score. Deleterious threshold > 0.5.",
				"##DESCRIPTION=MaxEntScan=MaxEntScan splicing prediction (reference bases score/alternate bases score).",
				"##DESCRIPTION=GeneSplicer=GeneSplicer splicing prediction (state/type/coordinates/confidence/score).",
				"##DESCRIPTION=dbscSNV=dbscSNV splicing prediction (ADA/RF score).",
				"##DESCRIPTION=COSMIC=COSMIC somatic variant database anntotation.",
				"##DESCRIPTION=NGSD_hom=Homozygous variant count in NGSD.",
				"##DESCRIPTION=NGSD_het=Heterozygous variant count in NGSD.",
				"##DESCRIPTION=NGSD_group=Homozygous / heterozygous variant count in NGSD with the same disease group (Neoplasms).",
				"##DESCRIPTION=classification=Classification from the NGSD.",
				"##DESCRIPTION=classification_comment=Classification comment from the NGSD.",
				"##DESCRIPTION=validation=Validation information from the NGSD. Validation results of other samples are listed in brackets!",
				"##DESCRIPTION=comment=Variant comments from the NGSD.",
				"##DESCRIPTION=gene_info=Gene information from NGSD (inheritance mode, gnomAD o/e scores).",
				"##FILTER=gene_blacklist=The gene(s) are contained on the blacklist of unreliable genes.",
				"##FILTER=off-target=Variant marked as 'off-target'.",
				"#chr	start	end	ref	obs	{$name}	filter	quality	gene	variant_type	coding_and_splicing	regulatory	OMIM	ClinVar	HGMD	RepeatMasker	dbSNP	1000g	gnomAD	gnomAD_hom_hemi	gnomAD_sub	phyloP	Sift	PolyPhen	fathmm-MKL	CADD	REVEL	MaxEntScan	GeneSplicer	dbscSNV	COSMIC	NGSD_hom	NGSD_het	NGSD_group	classification	classification_comment	validation	comment	gene_info",
			);
			file_put_contents($varfile, implode("\n", $content));
		}
		
	}

	// annotate CNV file
	if (file_exists($cnvfile))
	{
		$repository_basedir = repository_basedir();
		$data_folder = get_path("data_folder");
		if ($is_wes || $is_wgs || $is_wgs_shallow) //Genome/Exome: ClinCNV
		{

			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$cnvfile}", true);
			
		}
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$cnvfile}", true);
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes.bed -no_duplicates -url_decode -out {$cnvfile}", true);
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2021-10.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$cnvfile}", true);


		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2021_3.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$cnvfile}", true);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$cnvfile} -in2 {$omim_file} -no_duplicates -url_decode -out {$cnvfile}", true);
		}

		//annotate additional gene info
		$parser->exec(get_path("ngs-bits")."CnvGeneAnnotation", "-in {$cnvfile} -out {$cnvfile}", true);
		// skip annotation if no connection to the NGSD is possible
		if (db_is_enabled("NGSD"))
		{
			//annotate overlap with pathogenic CNVs
			$parser->exec(get_path("ngs-bits")."NGSDAnnotateCNV", "-in {$cnvfile} -out {$cnvfile}", true);
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
		//SV calling with manta
		$manta_evidence_dir = "{$folder}/manta_evid";
		create_directory($manta_evidence_dir);


		$manta_args = [
			"-bam ".$local_bamfile,
			"-evid_dir ".$manta_evidence_dir,
			"-out ".$sv_manta_file,
			"-threads ".$threads,
			"-build ".$sys['build'],
			"--log ".$parser->getLogFile()
		];
		if($has_roi) $manta_args[] = "-target ".$sys['target_file'];
		if(!$is_wgs) $manta_args[] = "-exome";
		
		$parser->execTool("NGS/vc_manta.php", implode(" ", $manta_args));

		// Rename Manta evidence file
		rename("$manta_evidence_dir/evidence_0.$name.bam", "$manta_evidence_dir/{$name}_manta_evidence.bam");
		rename("$manta_evidence_dir/evidence_0.$name.bam.bai", "$manta_evidence_dir/{$name}_manta_evidence.bam.bai");

		//create BEDPE files
		$parser->exec("{$ngsbits}VcfToBedpe", "-in $sv_manta_file -out $bedpe_out", true);
	}

	//add gene info annotation and NGSD counts
	if (db_is_enabled("NGSD"))
	{
		$parser->exec("{$ngsbits}BedpeGeneAnnotation", "-in $bedpe_out -out $bedpe_out -add_simple_gene_names", true);
		$parser->exec("{$ngsbits}NGSDAnnotateSV", "-in $bedpe_out -out $bedpe_out -ps $name", true);
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
}

//repeat expansions
if (in_array("sv", $steps) && !$annotation_only)
{
	//perform repeat expansion analysis (only for WGS/WES):
	$parser->execTool("NGS/vc_expansionhunter.php", "-in $local_bamfile -out $expansion_hunter_file -build ".$sys['build']." -pid $name");
}

//import to database
if (in_array("db", $steps))
{
	//import ancestry
	if(file_exists($ancestry_file)) $parser->execTool("NGS/db_import_ancestry.php", "-id {$name} -in {$ancestry_file}");
	
	//import QC
	$qc_files = array($qc_fastq, $qc_map);
	if (file_exists($qc_vc)) $qc_files[] = $qc_vc; 
	$parser->execTool("NGS/db_import_qc.php", "-id $name -files ".implode(" ", $qc_files)." -force");
	
	//check gender
	if(!$somatic) $parser->execTool("NGS/db_check_gender.php", "-in $bamfile -pid $name");	
	//import variants
	$args = ["-ps {$name}"];
	$import = false;
	if (file_exists($varfile) && !$is_wgs_shallow)
	{
		$args[] = "-var {$varfile}";
		$args[] = "-var_force";
		$import = true;
	}
	if (file_exists($cnvfile))
	{
		$args[] = "-cnv {$cnvfile}";
		$args[] = "-cnv_force";
		$import = true;
	}
	if (file_exists($bedpe_out))
	{
		$args[] = "-sv {$bedpe_out}";
		$args[] = "-sv_force";
		$import = true;
	}
	if ($import)
	{
		$parser->exec("{$ngsbits}NGSDAddVariantsGermline", implode(" ", $args), true);
	}
}


// Create Circos plot only if variant or copy-number calling was done
if (in_array("vc", $steps) || in_array("cn", $steps))
{
	if ($is_wes || $is_wgs || $is_wgs_shallow)
	{
		if (file_exists($cnvfile))
		{
			if (file_exists($cnvfile2))
			{
				$parser->execTool("NGS/create_circos_plot.php", "-folder $folder -name $name -build ".$sys['build']);
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

?>
