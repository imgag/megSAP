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
$steps_all = array("ma", "vc", "cn", "sv", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, vc=variant calling, cn=copy-number analysis, sv=structural-variant analysis, db=import into NGSD.", true, "ma,vc,cn,sv,db");
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("annotation_only", "Performs only a reannotation of the already created variant calls.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
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
$has_roi = $sys['target_file']!="";
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
// $local_bamfile = $parser->tempFolder("local_bam")."/".$name.".bam"; //local copy of BAM file to reduce IO over network
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
	// BAM file path to convert to FastQ
	$bamfile_to_convert = $bamfile;

	//determine input FASTQ files
	$fastq_regex = $folder."/*.fastq.gz";
	$fastq_files = glob($fastq_regex);
	

	if (count($fastq_files)==0)
	{
		if(file_exists($bamfile_to_convert))
		{
			trigger_error("No FASTQ files found in folder. Using BAM file to generate FASTQ files: $bamfile_to_convert", E_USER_NOTICE);

			trigger_error("Implement BamToFastq for Longreads!", E_USER_ERROR);
		}
		else
		{
			trigger_error("Found no read files found matching '$in_for' or '$in_rev'!", E_USER_ERROR);
		}
	}
	
	// run mapping
	$parser->execTool("NGS/mapping_minimap.php", "-in ".implode(" ", $fastq_files)." -out {$bamfile} -sample {$name} -threads {$threads} -system {$system}");

}
else
{
	//check genome build of BAM
	check_genome_build($bamfile, $build);
	
	$local_bamfile = $bamfile;
}

//variant calling
if (in_array("vc", $steps))
{
	// skip VC if only annotation should be done
	if (!$annotation_only)
	{		
		//TODO: implement
		$args = [];
		$args[] = "-bam ".$bamfile;
		$args[] = "-folder ".$folder;
		$args[] = "-name ".$name;
		$args[] = "-target ".$sys['target_file'];
		$args[] = "-target_extend 200";
		$args[] = "-threads ".$threads;
		$args[] = "-build ".$build;
		$args[] = "--log ".$parser->getLogFile();

		# determine model
		if($sys["name_short"] == "SQK-LSK114")
		{
			$args[] = "-model ".get_path("clair3_models")."/r1041_e82_400bps_hac_g632/";
		}
		else if($sys["name_short"] == "SQK-LSK109")
		{
			$args[] = "-model ".get_path("clair3_models")."/r941_prom_sup_g5014/";
		}
		else
		{
			trigger_error("Unsupported processing system '".$sys["shortname"]."' provided!", E_USER_ERROR);
		}

		
		$parser->execTool("NGS/vc_clair.php", implode(" ", $args));

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

	//TODO: test/implement
	$parser->execTool("Pipelines/annotate.php", implode(" ", $args));

	
	//determine ancestry
	if (ngsbits_build($build) != "non_human")
	{
		$parser->exec($ngsbits."SampleAncestry", "-in {$vcffile} -out {$ancestry_file} -build ".ngsbits_build($build), true);
	}
}

// //copy-number analysis
// if (in_array("cn", $steps))
// {
// 	// skip CN calling if only annotation should be done
// 	if (!$annotation_only)
// 	{
// 		//TODO: implement
// 	}
// 	else
// 	{
// 		check_genome_build($cnvfile, $build);
// 	}

// 	// annotate CNV file
// 	if (file_exists($cnvfile))
// 	{
// 		//TODO: implement
// 	}
// 	else
// 	{
// 		trigger_error("CNV file {$cnvfile} does not exist, skipping CNV annotation!", E_USER_WARNING);
// 	}
// }

// //structural variants
// if (in_array("sv", $steps))
// {
// 	// skip SV calling if only annotation should be done	
// 	if (!$annotation_only)
// 	{
		
		
// 		//create BEDPE files
// 		$parser->exec("{$ngsbits}VcfToBedpe", "-in $sv_manta_file -out $bedpe_out", true);
// 	}
// 	else
// 	{
// 		check_genome_build($sv_manta_file, $build);
// 	}


// 	//add gene info annotation
// 	if (db_is_enabled("NGSD"))
// 	{
// 		$parser->exec("{$ngsbits}BedpeGeneAnnotation", "-in $bedpe_out -out $bedpe_out -add_simple_gene_names", true);
// 	}

// 	//add NGSD counts from flat file
// 	$ngsd_annotation_folder = get_path("data_folder")."/dbs/NGSD/";
// 	$ngsd_sv_files = array("sv_deletion.bedpe.gz", "sv_duplication.bedpe.gz", "sv_insertion.bedpe.gz", "sv_inversion.bedpe.gz", "sv_translocation.bedpe.gz");
// 	$db_file_dates = array();

// 	// check file existance
// 	$all_files_available = file_exists($ngsd_annotation_folder."sv_breakpoint_density.igv");
// 	foreach ($ngsd_sv_files as $filename) 
// 	{
// 		if(!(file_exists($ngsd_annotation_folder.$filename) && file_exists($ngsd_annotation_folder.$filename.".tbi")))
// 		{
// 			$all_files_available = false;
// 			break;
// 		}
// 	}
// 	if ($all_files_available)
// 	{
// 		// store flat file modification date to detect changes during annotation 
// 		foreach ($ngsd_sv_files as $filename)
// 		{
// 			$db_file_dates[$filename] = filemtime($ngsd_annotation_folder.$filename);
// 			if ($db_file_dates[$filename] == false)
// 			{
// 				trigger_error("Cannot get modification date of '".$ngsd_annotation_folder.$filename."'!",E_USER_ERROR);
// 			}
// 		}
		
// 		//perform annotation
// 		$parser->exec("{$ngsbits}BedpeAnnotateCounts", "-in $bedpe_out -out $bedpe_out -processing_system ".$sys["name_short"]." -ann_folder {$ngsd_annotation_folder}", true);
// 		$parser->exec("{$ngsbits}BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv", true);

// 		// check if files changed during annotation
// 		foreach ($ngsd_sv_files as $filename)
// 		{
// 			if ($db_file_dates[$filename] != filemtime($ngsd_annotation_folder.$filename))
// 			{
// 				trigger_error("Annotation file '".$ngsd_annotation_folder.$filename."' has changed during annotation!",E_USER_ERROR);
// 			}
// 		}

// 	}
// 	else
// 	{
// 		trigger_error("Cannot annotate NGSD counts! At least one required file in '{$ngsd_annotation_folder}' is missing!", E_USER_WARNING);
// 	}
	

// 	//add optional OMIM annotation

// 	$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; 
// 	if(file_exists($omim_file))//OMIM annotation (optional because of license)
// 	{
// 		$parser->exec("{$ngsbits}BedpeAnnotateFromBed", "-in $bedpe_out -out $bedpe_out -bed $omim_file -url_decode -replace_underscore -col_name OMIM", true);
// 	}

// 	//add CNV overlap annotation
// 	if (file_exists($cnvfile))
// 	{
// 		$parser->exec("{$ngsbits}BedpeAnnotateCnvOverlap", "-in $bedpe_out -out $bedpe_out -cnv $cnvfile", true);
// 	}
// }

// // create Circos plot - if small variant, CNV or SV calling was done
// if ((in_array("vc", $steps) || in_array("cn", $steps) || in_array("sv", $steps)) && !$annotation_only)
// {
// 	if ($is_wes || $is_wgs || $is_wgs_shallow)
// 	{
// 		if (file_exists($cnvfile))
// 		{
// 			if (file_exists($cnvfile2))
// 			{
// 				$parser->execTool("NGS/create_circos_plot.php", "-folder $folder -name $name -build ".$build);
// 			}
// 			else
// 			{
// 				trigger_error("CNV file $cnvfile2 missing. Cannot create Circos plot!", E_USER_WARNING);
// 			}
// 		}
// 		else
// 		{
// 			trigger_error("CNV file $cnvfile missing. Cannot create Circos plot!", E_USER_WARNING);
// 		}
// 	}
// }

// // collect other QC terms - if CNV or SV calling was done
// if ((in_array("cn", $steps) || in_array("sv", $steps)) && !$annotation_only)
// {
// 	$terms = [];
// 	$sources = [];
	
// 	//CNVs
// 	if (file_exists($cnvfile))
// 	{
// 		$cnv_count_hq = 0;
// 		$cnv_count_hq_autosomes = 0;
// 		$cnv_count_loss = 0;
// 		$cnv_count_gain = 0;
// 		$h = fopen2($cnvfile, 'r');
// 		while(!feof($h))
// 		{
// 			$line = trim(fgets($h));
// 			if ($line=="") continue;
			
// 			if (starts_with($line, "##mean correlation to reference samples:"))
// 			{
// 				$value = trim(explode(":", $line)[1]);
// 				$terms[] = "QC:2000114\t{$value}";
// 			}
// 			if (starts_with($line, "##number of iterations:"))
// 			{
// 				$value = trim(explode(":", $line)[1]);
// 				$terms[] = "QC:2000115\t{$value}";
// 			}
			
// 			if ($line[0]!="#")
// 			{
// 				$parts = explode("\t", $line);
// 				$ll = $parts[4];
// 				if ($ll>=20)
// 				{
// 					++$cnv_count_hq;
					
// 					$chr = $parts[0];
// 					if (is_numeric(strtr($chr, ["chr"=>""])))
// 					{
// 						++$cnv_count_hq_autosomes;
// 						$cn = $parts[3];
// 						if ($cn<2) ++$cnv_count_loss;
// 						if ($cn>2) ++$cnv_count_gain;
// 					}
// 				}
// 			}
// 		}
// 		fclose($h);
		
// 		//counts (all, loss, gain)
// 		$terms[] = "QC:2000113\t{$cnv_count_hq}";
// 		if ($cnv_count_hq_autosomes>0)
// 		{
// 			$terms[] = "QC:2000118\t".number_format(100.0*$cnv_count_loss/$cnv_count_hq_autosomes, 2);
// 			$terms[] = "QC:2000119\t".number_format(100.0*$cnv_count_gain/$cnv_count_hq_autosomes, 2);
// 		}
		
// 		$sources[] = $cnvfile;
// 	}

// 	//SVs
// 	if (file_exists($bedpe_out))
// 	{
// 		$sv_count_pass = 0;
// 		$sv_count_del = 0;
// 		$sv_count_dup = 0;
// 		$sv_count_ins = 0;
// 		$sv_count_inv = 0;
// 		$sv_count_bnd = 0;
// 		$h = fopen2($bedpe_out, 'r');
// 		while(!feof($h))
// 		{
// 			$line = trim(fgets($h));
// 			if ($line=="" || $line[0]=="#") continue;
			
			
// 			$parts = explode("\t", $line);
// 			$filter = trim($parts[11]);
// 			if ($filter=="PASS")
// 			{
// 				++$sv_count_pass;
// 				$type = trim($parts[10]);
// 				if ($type=="DEL") ++$sv_count_del;
// 				if ($type=="DUP") ++$sv_count_dup;
// 				if ($type=="INS") ++$sv_count_ins;
// 				if ($type=="INV") ++$sv_count_inv;
// 				if ($type=="BND") ++$sv_count_bnd;
				
// 			}
// 		}
// 		fclose($h);
		
// 		$terms[] = "QC:2000117\t{$sv_count_pass}";
// 		if ($sv_count_pass>0)
// 		{
// 			$terms[] = "QC:2000120\t".number_format(100.0*$sv_count_del/$sv_count_pass, 2);
// 			$terms[] = "QC:2000121\t".number_format(100.0*$sv_count_dup/$sv_count_pass, 2);
// 			$terms[] = "QC:2000122\t".number_format(100.0*$sv_count_ins/$sv_count_pass, 2);
// 			$terms[] = "QC:2000123\t".number_format(100.0*$sv_count_inv/$sv_count_pass, 2);
// 			$terms[] = "QC:2000124\t".number_format(100.0*$sv_count_bnd/$sv_count_pass, 2);
// 		}
		
// 		$sources[] = $bedpe_out;
// 	}
	
// 	//create qcML file
// 	$tmp = $parser->tempFile("qc.tsv");
// 	file_put_contents($tmp, implode("\n", $terms));
// 	$parser->exec("{$ngsbits}TsvToQC", "-in $tmp -out $qc_other -sources ".implode(" ", $sources));
// }

// //import to database
// if (in_array("db", $steps))
// {
// 	//import ancestry
// 	if(file_exists($ancestry_file)) $parser->execTool("NGS/db_import_ancestry.php", "-id {$name} -in {$ancestry_file}");
	
// 	//import QC
// 	$qc_files = array($qc_fastq, $qc_map);
// 	if (file_exists($qc_vc)) $qc_files[] = $qc_vc; 
// 	if (file_exists($qc_other)) $qc_files[] = $qc_other; 
// 	$parser->execTool("NGS/db_import_qc.php", "-id $name -files ".implode(" ", $qc_files)." -force");
	
// 	//check gender
// 	if(!$somatic) $parser->execTool("NGS/db_check_gender.php", "-in $bamfile -pid $name");	
// 	//import variants
// 	$args = ["-ps {$name}"];
// 	$import = false;
// 	if (file_exists($varfile) && !$is_wgs_shallow)
// 	{
// 		//check genome build
// 		check_genome_build($varfile, $build);
		
// 		$args[] = "-var {$varfile}";
// 		$args[] = "-var_force";
// 		$import = true;
// 	}
// 	if (file_exists($cnvfile))
// 	{
// 		//check genome build
// 		//this is not possible for CNVs because the file does not contain any information about it
		
// 		$args[] = "-cnv {$cnvfile}";
// 		$args[] = "-cnv_force";
// 		$import = true;
// 	}
// 	if (file_exists($bedpe_out))
// 	{
// 		//check genome build
// 		check_genome_build($bedpe_out, $build);
		
// 		$args[] = "-sv {$bedpe_out}";
// 		$args[] = "-sv_force";
// 		$import = true;
// 	}
// 	if ($import)
// 	{
// 		$parser->exec("{$ngsbits}NGSDAddVariantsGermline", implode(" ", $args), true);
// 	}
// }

?>
