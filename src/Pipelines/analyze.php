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
$steps_all = array("ma", "vc", "an", "cn", "db", "sv");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, an=annotation, db=import into NGSD, cn=copy-number analysis, sv=structural-variant analysis.", true, "ma,vc,an,cn,db,sv");
$parser->addFlag("backup", "Backup old analysis files to old_[date] folder.");
$parser->addFlag("lofreq", "Add low frequency variant detection.", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("clip_overlap", "Soft-clip overlapping read pairs.", true);
$parser->addFlag("no_abra", "Skip realignment with ABRA.", true);
$parser->addFlag("correction_n", "Use Ns for errors by barcode correction.", true);
$parser->addString("out_folder", "Folder where analysis results should be stored. Default is same as in '-folder' (e.g. Sample_xyz/).", true, "default");
$parser->addFlag("somatic", "Set somatic single sample analysis options (i.e. correction_n, clip_overlap).");
extract($parser->parse($argv));

//init
if($out_folder=="default")
{
	$out_folder = $folder;
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

//set up local NGS data copy (to reduce network traffic and speed up analysis)
$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);

//output file names
//rename out folder
$bamfile = $out_folder."/".$name.".bam";
if(!in_array("ma", $steps))	$bamfile = $folder."/".$name.".bam";
$vcffile = $out_folder."/".$name."_var.vcf.gz";
$vcffile_annotated = $out_folder."/".$name."_var_annotated.vcf.gz";	
if(!in_array("vc", $steps))	$vcffile = $folder."/".$name."_var.vcf.gz";
$varfile = $out_folder."/".$name.".GSvar";
$lowcov_file = $out_folder."/".$name."_".$sys["name_short"]."_lowcov.bed";
if(!in_array("an", $steps))	$varfile = $folder."/".$name.".GSvar";
$log_ma  = $out_folder."/".$name."_log1_map.log";
$log_vc  = $out_folder."/".$name."_log2_vc.log";
$log_an  = $out_folder."/".$name."_log3_anno.log";
$log_db  = $out_folder."/".$name."_log4_db.log";
$log_cn  = $out_folder."/".$name."_log5_cn.log";
$log_sv = $out_folder ."/".$name."_log6_sv.log";
$qc_fastq  = $out_folder."/".$name."_stats_fastq.qcML";
$qc_map  = $out_folder."/".$name."_stats_map.qcML";
$qc_vc  = $out_folder."/".$name."_stats_vc.qcML";
if ($is_wgs || $is_wgs_shallow || $is_wes) //Genome/Exome: ClinCNV
{
	$cnvfile = $out_folder."/".$name."_cnvs_clincnv.tsv";
	$cnvfile2 = $out_folder."/".$name."_cnvs_clincnv.seg";
}
else //Panel: CNVHunter
{
	$cnvfile = $out_folder."/".$name."_cnvs.tsv";
	$cnvfile2 = $out_folder."/".$name."_cnvs.seg";
}	
$rohfile = $out_folder."/".$name."_rohs.tsv";
$baffile = $out_folder."/".$name."_bafs.igv";

$sv_manta_file = $out_folder ."/". $name . "_manta_var_structural.vcf.gz";
$small_indel_manta_file =  $out_folder ."/". $name . "_manta_var_smallIndels.vcf.gz";

//move old data to old_[date]_[random]-folder
if($backup && in_array("ma", $steps))
{
	$backup_pattern = "$name*.*,analyze*.log";
	$skip_pattern = array();
	$skip_pattern[] = "\w+\.fastq\.gz$";
	$skip_pattern[] = "SampleSheet\.csv$";
	if(!is_null($parser->getLogFile()))	$skip_pattern[] = $parser->getLogFile()."$";
	backup($out_folder, $backup_pattern, "#".implode("|", $skip_pattern)."#");
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
	if (count($files1)==0)
	{
		if(file_exists($bamfile))
		{
			trigger_error("No FASTQ files found in folder. Using BAM file to generate FASTQ files.", E_USER_NOTICE);

			// extract reads from BAM file
			$in_fq_for = $folder."/{$name}_BamToFastq_R1_001.fastq.gz";
			$in_fq_rev = $folder."/{$name}_BamToFastq_R2_001.fastq.gz";
			$parser->exec("{$ngsbits}BamToFastq", "-in $bamfile -out1 $in_fq_for -out2 $in_fq_rev -remove_duplicates", true);

			// use generated fastq files for mapping
			$files1 = array($in_fq_for);
			$files2 = array($in_fq_rev);
		}
		else
		{
			trigger_error("Found no read files found matching '$in_for' or '$in_rev'!", E_USER_ERROR);
		}
		
	}
	
	$args = array();
	if($clip_overlap) $args[] = "-clip_overlap";
	if($no_abra) $args[] = "-no_abra";
	if(file_exists($log_ma)) unlink($log_ma);
	if($correction_n) $args[] = "-correction_n";
	if(!empty($files_index)) $args[] = "-in_index " . implode(" ", $files_index);
	$parser->execTool("Pipelines/mapping.php", "-in_for ".implode(" ", $files1)." -in_rev ".implode(" ", $files2)." -system $system -out_folder $out_folder -out_name $name --log $log_ma ".implode(" ", $args)." -threads $threads");

	//low-coverage report
	if ($has_roi && !$is_wgs && !$is_wgs_shallow)
	{	
		$parser->exec("{$ngsbits}BedLowCoverage", "-in ".$sys['target_file']." -bam $bamfile -out $lowcov_file -cutoff 20", true);
		if (db_is_enabled("NGSD"))
		{
			$parser->exec("{$ngsbits}BedAnnotateGenes", "-in $lowcov_file -clear -extend 25 -out $lowcov_file", true);
		}
	}

	// delete fastq files after mapping
	$delete_fastq_files = get_path("delete_fastq_files", true);
	if ($delete_fastq_files==true || $delete_fastq_files=="true")
	{
		foreach($files1 as $fastq_file)
		{
			unlink($fastq_file);
		}
		foreach($files2 as $fastq_file)
		{
			unlink($fastq_file);
		}
	}
}

//variant calling
if (in_array("vc", $steps))
{
	if($is_wgs_shallow)
	{
		trigger_error("Skipping step 'vc' - Variant calling is not supported for shallow WGS samples!", E_USER_NOTICE);
	}
	else
	{
		$args = array();
		if ($has_roi)
		{
			$args[] = "-target ".$sys['target_file'];
			$args[] = "-target_extend 50";
		}
		if ($lofreq) //lofreq
		{
			$args[] = "-min_af 0.05";
		}
		else
		{
			$args[] = "-min_af 0.1";
		}
		if(file_exists($log_vc)) unlink($log_vc);
		
		//Do not call standard pipeline if there is only mitochondiral chrMT in target region
		$only_mito_in_target_region = false;
		if ($has_roi) $only_mito_in_target_region = exec2("cat ".$sys['target_file']." | cut -f1 | uniq")[0][0] == "chrMT";
		if(!$only_mito_in_target_region)
		{
			$parser->execTool("NGS/vc_freebayes.php", "-bam $bamfile -out $vcffile -build ".$sys['build']." --log $log_vc -threads $threads ".implode(" ", $args));
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
			$args[] = "-target $target_mito";
			$vcffile_mito = $parser->tempFile("_mito.vcf.gz");
			$parser->execTool("NGS/vc_freebayes.php", "-bam $bamfile -out $vcffile_mito -build ".$sys['build']." --log $log_vc ".implode(" ", $args));
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
		$params[] = "-bam {$bamfile}";
		$params[] = "-out {$baffile}";
		$params[] = "-build ".$sys['build'];
		if ($is_wgs)
		{
			$params[] = "-downsample 100";
		}
		$parser->execTool("NGS/baf_germline.php", implode(" ", $params));
	}
}

//annotation and reports
if (in_array("an", $steps))
{
	if($is_wgs_shallow)
	{
		trigger_error("Skipping step 'an' - Variant annotation is not supported for shallow WGS samples!", E_USER_NOTICE);
	}
	else
	{
		if(file_exists($log_an)) unlink($log_an);

		//annotate
		$args = array("-out_name $name", "-out_folder {$out_folder}", "-system {$system}", "--log {$log_an}");
		if ($is_wgs) $args[] = "-updown";
		$args[] = "-threads {$threads}";
		$parser->execTool("Pipelines/annotate.php", implode(" ", $args));
		
		//ROH detection
		if ($is_wes || $is_wgs)
		{
			$args = array();
			$args[] = "-in $vcffile_annotated";
			$args[] = "-out $rohfile";
			$args[] = "-var_af_keys AF,gnomAD_AF,gnomADg_AF"; //use 1000g, gnomAD exome, genomAD genome
			$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; //optional because of license
			$args[] = "-annotate ".repository_basedir()."/data/gene_lists/genes.bed ".(file_exists($omim_file) ? $omim_file : "");
			$parser->exec("{$ngsbits}RohHunter", implode(" ", $args), true);
		}
	}
}

//copy-number analysis
if (in_array("cn", $steps))
{
	if(!$has_roi)
	{
		trigger_error("Skipping step 'cn' - Copy number analysis is only supported for processing systems with target region BED file!", E_USER_NOTICE);
	}
	else
	{
		//remove log file
		if(file_exists($log_cn)) unlink($log_cn);
		
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
		
		if ($is_wes || $is_wgs || $is_wgs_shallow) //Genome/Exome: ClinCNV
		{
			//WGS: create folder for binned coverage data - if missing
			if ($is_wgs || $is_wgs_shallow)
			{
				$bin_size = 1000;
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
			
			//create BED file with GC and gene anntations - if missing
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
				$parser->exec("{$ngsbits}BedCoverage", "-min_mapq 0 -decimals 4 -bam {$bamfile} -in {$bed} -out {$cov_tmp}", true);
				
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
				"-cov_max ".($is_wgs ? "100" : "200"),
				"-max_cnvs ".($is_wgs ? "2000" : "200"),
				"--log {$log_cn}",
			);
			$parser->execTool("NGS/vc_clincnv_germline.php", implode(" ", $args), true);
			
			//copy results to output folder
			if (file_exists($cnv_out)) $parser->moveFile($cnv_out, $cnvfile);
			if (file_exists($cnv_out2)) $parser->moveFile($cnv_out2, $cnvfile2);
		}
		else //Panels: CnvHunter
		{
			//create coverage file
			$tmp_folder = $parser->tempFolder();
			$cov_file = $tmp_folder."/{$name}.cov";
			$parser->exec("{$ngsbits}BedCoverage", "-min_mapq 0 -decimals 4 -bam $bamfile -in ".$sys['target_file']." -out $cov_file", true);

			//copy coverage file to reference folder (has to be done before CnvHunter call to avoid analyzing the same sample twice)
			if (db_is_enabled("NGSD") && is_valid_ref_sample_for_cnv_analysis($name))
			{
				$ref_file = $ref_folder.$name.".cov";
				$parser->copyFile($cov_file, $ref_file);
				$cov_file = $ref_file;
			}
			
			//execute
			$cnv_out = $tmp_folder."/output.tsv";
			$cnv_out2 = $tmp_folder."/output.seg";
			$parser->execTool("NGS/vc_cnvhunter.php", "-min_z 3.5 -cov $cov_file -system $system -out $cnv_out -seg $name -n 20 --log $log_cn");

			//copy results to output folder
			if (file_exists($cnv_out)) $parser->moveFile($cnv_out, $cnvfile);
			if (file_exists($cnv_out2)) $parser->moveFile($cnv_out2, $cnvfile2);
		}
	}
}

//import to database
if (in_array("db", $steps))
{
	if(file_exists($log_db)) unlink($log_db);
	
	//import QC
	$qc_files = array($qc_fastq, $qc_map);
	if (file_exists($qc_vc)) $qc_files[] = $qc_vc; 
	$parser->execTool("NGS/db_import_qc.php", "-id $name -files ".implode(" ", $qc_files)." -force --log $log_db");
	
	//check gender
	$parser->execTool("NGS/db_check_gender.php", "-in $bamfile -pid $name --log $log_db");
	
	//import variants
	if (file_exists($varfile) || file_exists($cnvfile))
	{
		$args = ["-ps {$name}"];
		if (file_exists($varfile))
		{
			$args[] = "-var {$varfile}";
			$args[] = "-var_force";			
		}
		if (file_exists($cnvfile))
		{
			$args[] = "-cnv {$cnvfile}";
		}
		$parser->exec("{$ngsbits}NGSDAddVariantsGermline", implode(" ", $args), true);
	}
}

//structural variants
if (in_array("sv", $steps))
{
	if(!$is_wgs)
	{
		trigger_error("Skipping step 'sv' - Structural variant calling is only supported for WGS samples!", E_USER_NOTICE);
	}
	else
	{
		if(file_exists($log_sv)) unlink($log_sv);
		
		//SV calling with manta
		$manta_evidence_dir = "{$out_folder}/manta_evid";
		create_directory($manta_evidence_dir);

		$manta_args = [
			"-bam ".$bamfile,
			"-evid_dir ".$manta_evidence_dir,
			"-out ".$sv_manta_file,
			"-smallIndels ".$small_indel_manta_file,
			"-threads ".$threads,
			"-fix_bam",
			"-build ".$sys['build'],
			"--log ".$log_sv
		];
		if($has_roi) $manta_args[] = "-target ".$sys['target_file'];
		
		$parser->execTool("NGS/vc_manta.php", implode(" ", $manta_args));

		//create BEDPE files
		$bedpe_out = substr($sv_manta_file,0,-6)."bedpe";
		exec2("{$ngsbits}VcfToBedpe -in $sv_manta_file -out $bedpe_out");
	}
}

?>
