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
$parser->addInfile("system",  "Processing system INI file (determined from NGSD via the 'name' by default).", true);
$steps_all = array("ma", "vc", "an", "db", "cn");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, an=annotation, db=import into NGSD, cn=copy-number analysis.", true, implode(",", $steps_all));
$parser->addFlag("backup", "Backup old analysis files to old_[date] folder.");
$parser->addFlag("lofreq", "Add low frequency variant detection.", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addInt("thres", "Splicing region size used for annotation (flanking the exons).", true, 20);
$parser->addFlag("clip_overlap", "Soft-clip overlapping read pairs.", true);
$parser->addFlag("no_abra", "Skip realignment with ABRA.", true);
$parser->addString("out_folder", "Folder where analysis results should be stored. Default is same as in '-folder' (e.g. Sample_xyz/).", true, "default");
extract($parser->parse($argv));

$sys = load_system($system, $name);

if($out_folder=="default")
{
	$out_folder = $folder;
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
if(!in_array("vc", $steps))	$vcffile = $folder."/".$name."_var.vcf.gz";
$varfile = $out_folder."/".$name.".GSvar";
$lowcov_file = $out_folder."/".$name."_".$sys["name_short"]."_lowcov.bed";
if(!in_array("an", $steps))	$varfile = $folder."/".$name.".GSvar";
$log_ma  = $out_folder."/".$name."_log1_map.log";
$log_vc  = $out_folder."/".$name."_log2_vc.log";
$log_an  = $out_folder."/".$name."_log3_anno.log";
$log_db  = $out_folder."/".$name."_log4_db.log";
$log_cn  = $out_folder."/".$name."_log5_cn.log";
$qc_fastq  = $out_folder."/".$name."_stats_fastq.qcML";
$qc_map  = $out_folder."/".$name."_stats_map.qcML";
$qc_vc  = $out_folder."/".$name."_stats_vc.qcML";
$cnvfile = $out_folder."/".$name."_cnvs.tsv";
$cnvfile2 = $out_folder."/".$name."_cnvs.seg";

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
	$in_for = $folder."/*_R1_001.fastq.gz";
	$in_rev = $folder."/*_R2_001.fastq.gz";
	
	//find FastQ input files
	$files1 = glob($in_for);
	$files2 = glob($in_rev);
	if (count($files1)!=count($files2))
	{
		trigger_error("Found mismatching forward and reverse read file count!\n Forward: ".implode(" ", $in_for)."\n Reverse: ".implode(" ", $in_rev), E_USER_ERROR);
	}
	if (count($files1)==0)
	{
		trigger_error("Found no read files found matching '$in_for' or '$in_rev'!", E_USER_ERROR);
	}
	
	$args = array();
	if($clip_overlap) $args[] = "-clip_overlap";
	if($no_abra) $args[] = "-no_abra";
	if(file_exists($log_ma)) unlink($log_ma);
	
	$parser->execTool("Pipelines/mapping.php", "-in_for ".implode(" ", $files1)." -in_rev ".implode(" ", $files2)." -system $system -out_folder $out_folder -out_name $name --log $log_ma ".implode(" ", $args)." -threads $threads");
}

//variant calling
if (in_array("vc", $steps))
{	
	$args = array();
	if ($sys['target_file']!="")
	{
		$args[] = "-target ".$sys['target_file'];
		$args[] = "-target_extend 50";
	}
	if ($lofreq) //lofreq
	{
		$args[] = "-min_af 0.05";
	}
	else if (!$sys['shotgun']) //amplicon panels
	{
		$args[] = "-min_af 0.1";
	}
	
	if(file_exists($log_vc)) unlink($log_vc);
	$parser->execTool("NGS/vc_freebayes.php", "-bam $bamfile -out $vcffile -build ".$sys['build']." --log $log_vc ".implode(" ", $args));
	
	//if WES, perform special variant calling for mitochondria
	$mito = ($sys['type']=="WES");
	if ($mito)
	{
		$target_mito = $parser->tempFile("_mito.bed");
		if ($sys['build']=="hg19")
		{
			file_put_contents($target_mito, "chrM\t0\t16571");
		}
		else if ($sys['build']=="GRCh37")
		{
			file_put_contents($target_mito, "chrMT\t0\t16569");
		}
		else
		{
			trigger_error("No mitochondria target region available for genome ".$sys['build']."!", E_USER_ERROR);
		}
		$args = array();
		$args[] = "-no_ploidy";
		$args[] = "-min_af 0.01";
		$args[] = "-target $target_mito";
		$vcffile_mito = $parser->tempFile("_mito.vcf.gz");
		$parser->execTool("NGS/vc_freebayes.php", "-bam $bamfile -out $vcffile_mito -build ".$sys['build']." --log $log_vc ".implode(" ", $args));
	}
	
	//add sample header to VCF
	$hr = gzopen($vcffile, "r");
	if ($hr===FALSE) trigger_error("Could not open file '" + $vcffile + "'.", E_USER_ERROR);
	$vcf = $parser->tempFile("_unzipped.vcf");
	$hw = fopen($vcf, "w");
	if ($hw===FALSE) trigger_error("Could not open file '" + $vcf + "'.", E_USER_ERROR);
	while(!gzeof($hr))
	{
		$line = trim(gzgets($hr));
		if (strlen($line)==0) continue;
		if ($line[0]=="#" && $line[1]!="#")
		{
			fwrite($hw, "##ANALYSISTYPE=GERMLINE_SINGLESAMPLE\n");
			fwrite($hw, gsvar_sample_header($name, array("DiseaseStatus"=>"affected")));
		}
		fwrite($hw, $line."\n");
	}
	fclose($hr);
	if ($mito) //special variant calling for mitochondria 
	{
		$hr = gzopen($vcffile_mito, "r");
		if ($hr===FALSE) trigger_error("Could not open file '" + $vcffile_mito + "'.", E_USER_ERROR);
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
}

//annotation and reports
if (in_array("an", $steps))
{
	if(file_exists($log_an)) unlink($log_an);
	
	//annotate
	$args = array("-out_name $name", "-out_folder $out_folder", "-system $system", "-thres $thres", "--log $log_an");
	if (!db_is_enabled("NGSD")) $args[] = "-no_ngsd";
	$parser->execTool("Pipelines/annotate.php", implode(" ", $args));
	
	//low-coverage report
	if($sys['type']=="WGS") //WGS
	{
		$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-wgs -bam $bamfile -out $lowcov_file -cutoff 20", true);
		if (db_is_enabled("NGSD")) $parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $lowcov_file -extend 25 -out $lowcov_file", true);
	}
	else if ($sys['target_file']!="") //ROI (but not WGS)
	{	
		$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in ".$sys['target_file']." -bam $bamfile -out $lowcov_file -cutoff 20", true);
		if (db_is_enabled("NGSD")) $parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $lowcov_file -extend 25 -out $lowcov_file", true);
	}

	//x-diagnostics report
	if ($sys['name_short']=="ssX" || starts_with($sys['name_short'], "hpXLIDv"))
	{
		$parser->execTool("Pipelines/x_diagnostics.php", "-bam $bamfile -out_folder {$out_folder} --log $log_an -system {$system}");
	}
}

//import to database
if (in_array("db", $steps))
{
	if(file_exists($log_db)) unlink($log_db);
	$parser->execTool("NGS/db_check_gender.php", "-in $bamfile -pid $name --log $log_db");
	
	//import variants
	if (file_exists($varfile))
	{
		$parser->execTool("NGS/db_import_variants.php", "-id $name -var $varfile -build ".$sys['build']." -force --log $log_db");
	}
	
	//update last_analysis column of processed sample in NGSD (before db_import_qc.php because that can throw an error because of low coverage)
	updateLastAnalysisDate($name, $bamfile);

	//import QC data
	$qc_files = array($qc_fastq, $qc_map);
	if (file_exists($qc_vc)) $qc_files[] = $qc_vc; 
	$parser->execTool("NGS/db_import_qc.php", "-id $name -files ".implode(" ", $qc_files)." -force --log $log_db");
}

//copy-number analysis
if (in_array("cn", $steps) && $sys['type']!="WGS")
{
	if(file_exists($log_cn)) unlink($log_cn);
	
	//create coverage file
	$tmp_folder = $parser->tempFolder();
	$cov_file = $tmp_folder."/{$name}.cov";
	$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -bam $bamfile -in ".$sys['target_file']." -out $cov_file", true);

	//copy coverage file to reference folder (has to be done before CnvHunter call to avoid analyzing the same sample twice)
	if (is_valid_ref_sample_for_cnv_analysis($name))
	{
		//create reference folder if it does not exist
		$ref_folder = get_path("data_folder")."/coverage/".$sys['name_short']."/";
		if (!is_dir($ref_folder)) mkdir($ref_folder);
		
		//copy file
		$ref_file = $ref_folder.$name.".cov";
		copy2($cov_file, $ref_file);
		$cov_file = $ref_file;
	}
	
	//perform copy-number analysis
	$cnv_out = $tmp_folder."/output.tsv";
	$cnv_out2 = $tmp_folder."/output.seg";
	$parser->execTool("NGS/vc_cnvhunter.php", "-min_z 3.5 -cov $cov_file -system $system -out $cnv_out -seg $name --log $log_cn");

	//copy results to output folder
	if (file_exists($cnv_out)) copy2($cnv_out, $cnvfile);
	if (file_exists($cnv_out2)) copy2($cnv_out2, $cnvfile2);
}

