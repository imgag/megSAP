<?php

/**
	@page trio
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function determine_index($name, $parts)
{
	$index = array_search($name, $parts);
	
	if ($index===FALSE || $index==-1)
	{
		trigger_error("Could not determine index of column '$name' in header line: ".implode(" ", $parts), E_USER_ERROR);
	}
	
	return $index;
}

function fix_gsvar_file($gsvar, $sample_c, $sample_f, $sample_m, $gender_data)
{
	global $parser;
	
	$h = fopen2($gsvar, "r");
	$tmp = $parser->tempFile(".GSvar");
	$h2 = fopen2($tmp, "w");
	while(!feof($h))
	{
		$line = nl_trim(fgets($h));
		if ($line=="") continue;
		
		//compare child gender from data with gender from header
		if (starts_with($line, "##SAMPLE=<ID={$sample_c},"))
		{
			//determine gender from header
			$gender_header = "n/a";
			$parts = explode(",", substr($line, 10, -1));
			foreach($parts as $part)
			{
				if (starts_with($part, "Gender="))
				{
					list(, $gender_header) = explode("=", $part);
				}
			}
			print "Gender of child (from header): {$gender_header}\n";
			
			//deviating => error
			if ($gender_data!=$gender_header && $gender_header!="n/a" && $gender_header!="" && $gender_data!="n/a")
			{
				trigger_error("Gender of child from sample header '{$gender_header}' deviates from gender from data '{$gender_data}'", E_USER_ERROR);
			}
				
			//replace gender in header by gender from data
			if ($gender_data!=$gender_header)
			{
				$line = strtr($line, array("Gender={$gender_header}"=>"Gender={$gender_data}"));
			}
		}
		
		//add gender of father/mother if missing (i.e. when NGSD is not enabled)
		if (starts_with($line, "##SAMPLE=<ID={$sample_f},") && !contains($line, "Gender="))
		{
			$line = substr($line, 0, -1).",Gender=male>";
		}
		if (starts_with($line, "##SAMPLE=<ID={$sample_m},") && !contains($line, "Gender="))
		{
			$line = substr($line, 0, -1).",Gender=female>";
		}
		
		//update analysis type
		if (starts_with($line, "##ANALYSISTYPE="))
		{
			$line = "##ANALYSISTYPE=GERMLINE_TRIO";
		}
		
		fwrite($h2, "$line\n");
	}
	fclose($h2);
	fclose($h);
	$parser->moveFile($tmp, $gsvar);
}

//parse command line arguments
$parser = new ToolBase("trio", "Trio analysis pipeline.");
$parser->addInfile("f", "BAM file of father.", false, true);
$parser->addInfile("m", "BAM file of mother.", false, true);
$parser->addInfile("c", "BAM file of child (index).", false, true);
$parser->addInfile("out_folder", "Output folder name.", false);
//optional
$parser->addInfile("system",  "Processing system INI file used for all samples (automatically determined from NGSD if the basename of 'c' is a valid processed sample name).", true);
$steps_all = array("vc", "cn", "sv", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nvc=variant calling, cn=copy-number analysis, cn=copy-number analysis, sv=structural variant calling, db=database import.", true, implode(",", $steps_all));
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("no_check", "Skip gender check of parents and parent-child correlation check (otherwise done before variant calling)");
$parser->addFlag("annotation_only", "Performs only a reannotation of the already created variant calls.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");

extract($parser->parse($argv));

//init
$sample_c = basename2($c);
$sample_f = basename2($f);
$sample_m = basename2($m);

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($out_folder."/trio_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//file names
$gsvar = "{$out_folder}/trio.GSvar";
$vcf_all = "{$out_folder}/all.vcf.gz";
$cnv_multi = "{$out_folder}/trio_cnvs_clincnv.tsv";
$bedpe_out = "{$out_folder}/trio_var_structural_variants.bedpe";

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

if ($annotation_only)
{
	// check if required VC files for annotation is available and print warning otherwise
	if (in_array("vc", $steps) && !file_exists($vcf_all))
	{
		trigger_error("VCF for reannotation is missing. Skipping 'vc' step!", E_USER_WARNING);
		if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
	}

	if (in_array("cn", $steps) && !file_exists($cnv_multi))
	{
		$cnv_multi = "{$out_folder}/trio_cnvs.tsv";
		if (!file_exists($cnv_multi))
		{
			trigger_error("CN file for reannotation is missing. Skipping 'cn' step!", E_USER_WARNING);
			if (($key = array_search("cn", $steps)) !== false) unset($steps[$key]);
		}
	}

	if (in_array("sv", $steps) && !file_exists($bedpe_out))
	{
		trigger_error("BEDPE file for reannotation is missing. Skipping 'sv' step!", E_USER_WARNING);
		if (($key = array_search("sv", $steps)) !== false) unset($steps[$key]);
	}
}

//prepare multi-sample paramters
$sys = load_system($system, $sample_c); //required in case the the system is unset
$genome = genome_fasta($sys['build']);
$args_multisample = [
	"-bams $c $f $m",
	"-status affected control control",
	"-out_folder $out_folder",
	"-system $system",
	"-prefix trio",
	"-no_sync", //already done if needed
	"-threads $threads"
	];
	
if ($annotation_only) $args_multisample[] = "-annotation_only";

//check steps
$is_wgs_shallow = $sys['type']=="WGS (shallow)";
if ($is_wgs_shallow)
{
	if (in_array("vc", $steps))
	{
		trigger_error("Skipping step 'vc' - Variant calling is not supported for shallow WGS samples!", E_USER_NOTICE);
		if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
	}
}

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);
}

//pre-analysis checks
if (!$no_check && !$is_wgs_shallow && !$annotation_only)
{
	//check parent-child correlation
	$min_corr = 0.45;
	$output = $parser->execApptainer("ngs-bits", "SampleSimilarity", "-in $f $c -mode bam -ref {$genome} -build ".ngsbits_build($sys['build']), [$f, $c, $genome]);
	$correlation = explode("\t", $output[0][1])[3];
	if ($correlation<$min_corr)
	{
		trigger_error("The genotype correlation of father and child is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
	}
	$output = $parser->execApptainer("ngs-bits", "SampleSimilarity", "-in $m $c -mode bam -ref {$genome} -build ".ngsbits_build($sys['build']), [$m, $c, $genome]);
	$correlation = explode("\t", $output[0][1])[3];
	if ($correlation<$min_corr)
	{
		trigger_error("The genotype correlation of mother and child is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
	}
	
	//check gender of parents
	$parser->execTool("Tools/db_check_gender.php", " -in $f -pid $sample_f -gender male");
	$parser->execTool("Tools/db_check_gender.php", " -in $m -pid $sample_m -gender female");
}

//variant calling (and annotation)
if (in_array("vc", $steps))
{
	//variant calling with multi-sample pipeline
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps vc", true); 

	//determine mendelian error rate
	list($stdout) = $parser->execApptainer("ngs-bits", "TrioMendelianErrors", "-vcf {$vcf_all} -c {$sample_c} -f {$sample_f} -m {$sample_m}", [$vcf_all]);
	foreach($stdout as $line)
	{
		if (starts_with($line, "Mendelian error rate:"))
		{
			print trim($line)."\n";
		}
	}
	
	//determine gender of child
	list($stdout, $stderr) = $parser->execApptainer("ngs-bits", "SampleGender", "-method hetx -in $c -build ".ngsbits_build($sys['build'])." -ref {$genome}", [$c, $genome]);
	$gender_data = explode("\t", $stdout[1])[1];
	if ($gender_data!="male" && $gender_data!="female") $gender_data = "n/a";
	print "Gender of child (from data): {$gender_data}\n";
	
	//write output file
	fix_gsvar_file($gsvar, $sample_c, $sample_f, $sample_m, $gender_data);
	
	//double check modified output file
	$parser->execTool("Tools/check_tsv.php", "-in $gsvar -build ".$sys['build']);
}

//copy-number and UPD
if (in_array("cn", $steps))
{
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps cn", true);
	
	if (!$annotation_only)
	{
		//UPD detection
		if (!$is_wgs_shallow)
		{
			$args_upd = [
				"-in {$out_folder}/all.vcf.gz",
				"-c {$sample_c}",
				"-f {$sample_f}",
				"-m {$sample_m}",
				"-out {$out_folder}/trio_upd.tsv",
				];
			$in_files = [$out_folder, $c];
			$base = dirname($c)."/".basename2($c);
			if (file_exists("{$base}_cnvs.tsv"))
			{
				$args_upd[] = "-exclude {$base}_cnvs.tsv";
			}
			else if (file_exists("{$base}_cnvs_clincnv.tsv"))
			{
				$args_upd[] = "-exclude {$base}_cnvs_clincnv.tsv";
			}
			else
			{
				trigger_error("Child CNV file not found for UPD detection!", E_USER_WARNING);
			}
			$parser->execApptainer("ngs-bits", "UpdHunter", implode(" ", $args_upd), $in_files);
		}
		
		//update GSvar file (because 'an' is skipped)
		if ($is_wgs_shallow)
		{
			//determine gender of child
			list($stdout, $stderr) = $parser->execApptainer("ngs-bits", "SampleGender", "-method xy -in $c -build ".ngsbits_build($sys['build'])." -ref {$genome}", [$c, $genome]);
			$gender_data = explode("\t", $stdout[1])[1];
			if ($gender_data!="male" && $gender_data!="female") $gender_data = "n/a";
			print "Gender of child (from data): {$gender_data}\n";
		
			fix_gsvar_file($gsvar, $sample_c, $sample_f, $sample_m, $gender_data);

			//contamination check (only necessary for shallow - for deep genomes we have the check in VariantQC)
			if($lines = file($cnv_multi))
			{
				list($stdout, $stderr) = $parser->execApptainer("ngs-bits", "TrioMaternalContamination", "-bam_m $m -bam_f $f -bam_c $c -build ".ngsbits_build($sys['build']), [$m, $f, $c]);
				$trio_info = "##TrioMaternalContamination ".implode(" | ", $stdout)."\n";
				foreach($lines as $line_num => $line)
				{
					if(!starts_with($line, "##"))
					{
						array_splice($lines, $line_num, 0, $trio_info);
						break;
					}
				}
				file_put_contents($cnv_multi, $lines);
			}
		}
	}
}

//sv calling
if (in_array("sv", $steps))
{
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps sv", true);
}

//everything worked > import/update sample data in NGSD
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD", false);
	
	//gender father
	$info_f = get_processed_sample_info($db_conn, $sample_f);
	$s_id_f = $info_f['s_id'];
	if ($info_f['gender']!="male")
	$db_conn->executeStmt("UPDATE sample SET gender='male' WHERE id='{$s_id_f}'");
	
	//gender mother
	$info_m = get_processed_sample_info($db_conn, $sample_m);
	$s_id_m = $info_m['s_id'];
	$db_conn->executeStmt("UPDATE sample SET gender='female' WHERE id='{$s_id_m}'");
	
	//sample relations
	$info_c = get_processed_sample_info($db_conn, $sample_c);
	$s_id_c = $info_c['s_id'];
	$db_conn->executeStmt("INSERT IGNORE INTO `sample_relations`(`sample1_id`, `relation`, `sample2_id`) VALUES ({$s_id_f},'parent-child',{$s_id_c})");
	$db_conn->executeStmt("INSERT IGNORE INTO `sample_relations`(`sample1_id`, `relation`, `sample2_id`) VALUES ({$s_id_m},'parent-child',{$s_id_c})");
	
	//add secondary analysis (if missing)
	$parser->execTool("Tools/db_import_secondary_analysis.php", "-type 'trio' -gsvar {$gsvar}");
}

?>
