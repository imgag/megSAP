<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

/**
	@page init_cnv_ref_folder
*/

// parse command line arguments
$parser = new ToolBase("init_cnv_ref_folder", "Initialize coverage folder of a processing system with data from all valid samples.");
$parser->addString("name", "Processing system short name.", false);
$parser->addFlag("clear", "Remove existing coverage files.");
$parser->addFlag("tumor_only", "Keep only tumor samples, normally tumor samples are removed.");
$parser->addFlag("somatic","Create coverage files for tumor-normal samples, including off-target coverage files.");
$parser->addString("build","reference genome",true,"GRCh37");
extract($parser->parse($argv));

//check system
$db = DB::getInstance("NGSD");
$res = $db->executeQuery("SELECT * FROM processing_system WHERE name_short=:sys", array("sys"=>$name));
if (count($res)!=1)
{
	trigger_error("Invalid processing system short name '$name'! Valid names are:\n".implode("\n", $db->getValues("SELECT name_short FROM processing_system ORDER BY name_short ASC")), E_USER_ERROR);	
}
$roi = trim($res[0]['target_file']);
if ($roi=="")
{
	trigger_error("Processing system '$name' has no target region file annotated!", E_USER_ERROR);	
}
$type = $res[0]['type'];
if ($type=="WGS" || $type=="WGS (shallow)")
{
	trigger_error("This script does not support WGS processing systems!", E_USER_ERROR);	
}

//Make list of samples that have same processing system
$ngsbits = get_path("ngs-bits");
$sample_table_file = $parser->tempFile(".tsv");
$pipeline= [
	["{$ngsbits}NGSDExportSamples", "-system $name -add_path"],
	["{$ngsbits}TsvFilter", "-filter 'system_name_short is $name'"], //this is necessary because NGSDExportSamples performs fuzzy match
	["{$ngsbits}TsvSlice", "-cols name,path -out $sample_table_file"],
];
$parser->execPipeline($pipeline, "NGSD sample extraction", true);
$samples = file($sample_table_file);

if(!$somatic)
{
	//create/clear folder
	$ref_folder = get_path("data_folder")."/coverage/$name".($tumor_only ? "-tumor" : "")."/";		
	if (!is_dir($ref_folder))
	{
		mkdir($ref_folder);
		if (!chmod($ref_folder, 0777))
		{
			trigger_error("Could not change privileges of folder '{$ref_folder}'!", E_USER_ERROR);
		}
	}
	if ($clear)
	{
		exec2("rm -rf $ref_folder/*.cov");
	}

	//print output
	print "Processing system: $name\n";
	print "Target file: $roi\n";
	print "Coverage folder: $ref_folder\n";
	print "\n";

	//check samples
	
	$bams = array();
	$valid = 0;
	foreach($samples as $line)
	{
		$line = trim($line);
		if ($line=="" || $line=="#") continue;
		
		list($sample, $path) = explode("\t", $line);
		
		//check sample is valid
		if(!is_valid_ref_sample_for_cnv_analysis($sample, $tumor_only)) continue;
		++$valid;
		
		//check bam
		$bam = "{$path}/{$sample}.bam";
		if (file_exists($bam))
		{
			$bams[] = $bam;
		}
	}
	print "Found ".count($samples)." samples for this processing system in NGSD.\n";
	print "Found $valid samples that are valid reference samples.\n";
	print "Found ".count($bams)." samples with BAM files.\n";

	//create new coverage files
	foreach($bams as $bam)
	{
		$cov_file = "$ref_folder/".basename($bam, ".bam").".cov";
		
		//skip existing coverage files
		if (file_exists($cov_file))
		{
			print "$bam skipped (coverage file already exists)\n";
			continue;
		}
		
		//calculate coverage file
		print "$bam processing...\n";
		exec2($ngsbits."BedCoverage -min_mapq 0 -decimals 4 -bam $bam -in $roi -out $cov_file");
	}

	//chmod
	exec2("chmod 775 $ref_folder/*.cov");
}
else //somatic tumor-normal pairs
{
	//normal coverage output files
	$off_target_bed = get_path("data_folder")."/coverage/off_target_beds/{$name}.bed";
	$ref_n_dir = get_path("data_folder")."/coverage/{$name}";
	if(!is_dir($ref_n_dir)) create_directory($ref_n_dir);
	$ref_off_n_dir = get_path("data_folder")."/coverage/{$name}_off_target";
	if(!is_dir($ref_off_n_dir)) create_directory($ref_off_n_dir);

	//tumor coverage output files
	$ref_t_dir = get_path("data_folder")."/coverage/{$name}-tumor";
	if(!is_dir($ref_t_dir)) create_directory($ref_t_dir);
	$ref_off_t_dir = get_path("data_folder")."/coverage/{$name}-tumor_off_target";
	if(!is_dir($ref_off_t_dir)) create_directory($ref_off_t_dir);
	
	
	print "Processing system: $name\n";
	print "Target file: $roi\n";
	print "Normal Coverage folder: $ref_n_dir\n";
	print "Normal off-target coverage folder: $ref_off_n_dir\n";
	print "Tumor Coverage folder: $ref_t_dir\n";
	print "Tumor off-target coverage folder: $ref_off_t_dir\n";
	print "\n";
	
	//list containing ids of corresponding tumor normal ids
	$t_n_list = "{$ref_t_dir}/list_tid-nid.csv";
	
	if($clear)
	{
		print "Clearing coverage file directories...\n\n";
		if(file_exists($off_target_bed)) exec2("rm  $off_target_bed");
		exec2("rm -rf {$ref_n_dir}/*.cov");
		exec2("rm -rf {$ref_off_n_dir}/*.cov");
		exec2("rm -rf {$ref_t_dir}/*.cov");
		if(file_exists($t_n_list))exec2("rm $t_n_list");
		exec2("rm -rf {$ref_off_t_dir}/*.cov");
	}
	
	if(!file_exists($off_target_bed))
	{
		print "Creating off-target BED-file {$off_target_bed}\n\n";
		create_off_target_bed_file($off_target_bed,$roi,genome_fasta($build));
	}
	
	
	$t_count = 0;
	$t_off_count = 0;
	$n_count = 0;
	$n_off_count = 0;
	
	foreach($samples as $line)
	{
		if(starts_with($line,"#")) continue;
		list($sample, $path) = explode("\t", trim($line));
		$info = get_processed_sample_info($db,$sample);

		$bam = $info["ps_bam"];
		
		if(!file_exists($bam))
		{
			print "Skipping $sample because $bam does not exist.\n";
			continue;
		}

		if($info['is_tumor'] == '1')
		{
			echo "tumor $bam processing...\n";
			if(!is_valid_ref_tumor_sample_for_cnv_analysis($sample))  continue;
			if(!file_exists("{$ref_t_dir}/{$sample}.cov"))
			{
				exec2($ngsbits."BedCoverage -min_mapq 0 -decimals 4 -bam $bam -in $roi -out {$ref_t_dir}/{$sample}.cov",true);
				++$t_count;
			}
			//off-target
			if(!file_exists("{$ref_off_t_dir}/{$sample}.cov"))
			{
				exec2($ngsbits."BedCoverage -min_mapq 10 -decimals 4 -in $off_target_bed -bam $bam -out {$ref_off_t_dir}/{$sample}.cov" ,true);
				++$t_off_count;
			}
		}
		else
		{
			echo "normal $bam processing...\n";
			if(!is_valid_ref_sample_for_cnv_analysis($sample)) continue;
			if(!file_exists("{$ref_n_dir}/{$sample}.cov"))
			{
				exec2($ngsbits."BedCoverage -min_mapq 0 -decimals 4 -bam $bam -in $roi -out {$ref_n_dir}/{$sample}.cov",true);
				++$n_count;
			}
			//off-target
			if(!file_exists("{$ref_off_n_dir}/{$sample}.cov"))
			{
				exec2($ngsbits."BedCoverage -min_mapq 10 -decimals 4 -in $off_target_bed -bam $bam -out {$ref_off_n_dir}/{$sample}.cov" ,true);
				++$n_off_count;
			}
		}
	}	
	
	print "\n";
	print "Processed $n_count normal samples and $t_count tumor samples\n";
	print "Processed $n_off_count normal off-target samples and $t_off_count tumor off-target samples\n";
	print "\n";
	
	
	//Create list with tumor normal pairs, has to be done every time any cov file is changed
	print "Creating list with tumor-normal ids {$t_n_list}\n";
	$handle = fopen2($t_n_list,"w");
	fputs($handle,"##THIS FILE CONTAINS TUMOR AND NORMAL IDS OF PROCESSING SYSTEM {$name}\n#tumor_id,normal_id\n");
	
	$t_covs = glob("{$ref_t_dir}/*.cov");
	
	foreach($t_covs as $t_cov)
	{
		$tid = substr(basename($t_cov),0,-4);
		$info = get_processed_sample_info($db,$tid,false);
		if($info === null) continue;
		$nid = $info['normal_name'];
		if(empty($nid) || !file_exists("{$ref_n_dir}/{$nid}.cov")) continue;
		fputs($handle,"{$tid},{$nid}\n");
	}
	
	fclose($handle);
	
}

?>
