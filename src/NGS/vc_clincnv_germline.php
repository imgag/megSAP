<?php 
/** 
	@page vc_clincnv_germline
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_clincnv_germline", "Germline copy-number calling with ClinCNV.");
$parser->addInfile("cov", "Coverage file for sample (tab-separated file with columns chr, start, end, coverage).", false);
$parser->addInfile("cov_folder", "Coverage files folder.", false);
$parser->addInfile("bed", "BED file with annotations e.g. GC-content and gene names.", false);
$parser->addOutFile("out", "Output file in TSV format.", false);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
//optional
$parser->addInt("cov_min", "Minimum number of coverage files required for CNV analysis.", true, 20);
$parser->addInt("cov_max", "Maximum number of coverage files used for CNV analysis, used to keep run-time and RAM requirement manageable (~200 for WES and ~100 for WGS).", true, 200);
$parser->addInt("max_cnvs", "Number of expected CNVs (~200 for WES and ~2000 for WGS).", true, 2000);
$parser->addInt("max_tries", "Maximum number of tries for calling ClinCNV (R parallelization sometimes breaks with no reason", true, 10);
$parser->addInt("regions", "Number of subsequent regions that must show a signal for a call.", true, 2);
$parser->addFlag("skip_super_recall", "Skip super-recall (down to one region and log-likelihood 3).");
//tumor only flags (use for somatic with no normal sample)
$parser->addFlag("tumor_only", "Analyze tumor sample without a paired normal sample.");
$parser->addInfile("bed_off","Off-target bed file.",true); //s_dna
$parser->addInfile("cov_off","Off-target coverage file for normal sample",true);
$parser->addString("cov_folder_off", "Folder with off-target normal coverage files (if different from [data_folder]/coverage/[system_short_name]_off_target/.", true, "auto");
$parser->addString("baf_folder","Folder containing files with B-Allele frequencies.",true);
$parser->addString("cov_off_min","Minimum number of off-target coverage files required for CNV analysis.",true, 10);

extract($parser->parse($argv));

//Creates file with paths to coverage files using ref folder and current sample cov path, if $sample_ids is set: skip all ids not contained
function create_file_with_paths($ref_cov_folder,$cov_path, &$sample_ids)
{
	global $parser;
	$new_sample_ids;

	$ref_paths = glob("{$ref_cov_folder}/*.cov");

	$paths_to_be_included = array();
	//Remove samples that are not in $sample_ids
	for($i=0;$i<count($ref_paths);++$i)
	{	
		$id = basename($ref_paths[$i], ".cov");

		if(in_array($id, $sample_ids))
		{
			$new_sample_ids[] = $id;
			$paths_to_be_included[] = $ref_paths[$i];
		}
	}

	//Check whether cov file is already in cov folder -> remove from list (could be specified in another dir!)
	for($i=0;$i<count($paths_to_be_included);++$i)
	{
		if(strpos($paths_to_be_included[$i],basename($cov_path,".cov")) !== false)
		{
			//for tumor only delete the old tumor cov file from folder if it exists, 
			//and add the new calculated cov file again and reassign correct indices
			unset($paths_to_be_included[$i]);
			$paths_to_be_included = array_values($paths_to_be_included);
			unset($new_sample_ids[$i]);
			$new_sample_ids = array_values($new_sample_ids);
			break;
		}
	}
	$paths_to_be_included[] = $cov_path;
	$new_sample_ids[] = basename($cov_path, "_off_target.cov");

	$out_file = $parser->tempFile(".txt");
	sort($paths_to_be_included);
	file_put_contents($out_file,implode("\n",$paths_to_be_included) );

	$sample_ids = $new_sample_ids;
	
	return $out_file;
}

function determine_rows_to_use($cov, $roi_nonpoly)
{	
	//determine non-polymorphic regions
	$nonpoly_regs = array();
	$file = file($roi_nonpoly);
	for($i=0; $i<count($file); ++$i)
	{
		$line = trim($file[$i]);
		if ($line=="" || $line[0]=="#") continue;
		
		list($chr, $start, $end) = explode("\t", $line);
		
		$nonpoly_regs["$chr:$start-$end"] = true;
	}
	
	//build bit array that determines which indices to use
	$output = array();	
	$file = file($cov);
	for($i=0; $i<count($file); ++$i)
	{
		$line = trim($file[$i]);
		
		//skip: empty/comment line
		if ($line=="" || $line[0]=="#")
		{
			$output[] = false;
			continue;
		}
		
		list($chr, $start, $end, $cov_avg) = explode("\t", $line);
		
		//skip: reference sample has no coverage
		if ($cov_avg==0.0) 
		{
			$output[] = false;
			continue;
		}
		
		//skip: polymorphic regions
		if (!isset($nonpoly_regs["$chr:$start-$end"]))
		{
			$output[] = false;
			continue;
		}
		
		//skip: all chromosomes but autosomes
		if (starts_with($chr, "chr")) $chr = substr($chr, 3);
		if (!is_numeric($chr))
		{
			$output[] = false;
			continue;
		}
				
		$output[] = true;
	}
	
	$output[] = false; //last emtpy line
	
	return $output;
}

//function to load coverage profiles
function load_coverage_profile($filename, &$rows_to_use, &$output)
{
	$output = array();

	$i = -1;
	$h = fopen2($filename, "r");
	while(!feof($h))
	{
		$line = fgets($h);
		
		++$i;
		if ($rows_to_use[$i]==false) continue;
		
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		$parts = explode("\t", $line);
		if (count($parts)<4) trigger_error("Invalid coverage file line in file {$filename}:\n{$line}", E_USER_ERROR);
		
		list($chr, $start, $end, $avg_cov) = $parts;
		$output[$chr][] = $avg_cov;
	}
}

//function to write an empty cnv file
function generate_empty_cnv_file($out, $command, $stdout, $ps_name, $error_message, $tumor_only)
{
	//ClinVAR did not generate CNV file
	//generate file with basic header lines
	$cnv_output = fopen($out, "w");
	if($tumor_only)
	{
		fwrite($cnv_output, "##ANALYSISTYPE=CLINCNV_TUMOR_ONLY\n");
	}
	else
	{
		fwrite($cnv_output, "##ANALYSISTYPE=CLINCNV_GERMLINE_SINGLE\n");
	}
	preg_match('/^.*ClinCNV-([\d\.]*).*$/', $command, $matches);
	if(sizeof($matches) >= 2)
	{
		fwrite($cnv_output, "##ClinCNV version: v{$matches[1]}\n");
	}
	fwrite($cnv_output, "##CNV calling skipped: {$error_message}.\n");
	//search in CLinCNV stdout for lines indicating failed QC
	foreach($stdout as $line)
	{
		//if QC for input sample failed, add it to header line
		if(strpos($line, "{$ps_name} did not pass QC") == true)
		{
			preg_match('/^.*\"(.*)\".*$/', $line, $matches);
			if(sizeof($mathces >= 2))
			{
				fwrite($cnv_output, "##{$matches[1]}\n");
			}
		}
	}
	fwrite($cnv_output, "#chr\tstart\tend\tCN_change\tloglikelihood\tno_of_regions\tlength_KB\tpotential_AF\tgenes\n");

	fclose($cnv_output);
}

//init
$repository_basedir = repository_basedir();
$ps_name = basename($cov,".cov");
$command = get_path("clincnv")."/clinCNV.R";

//determine coverage files
$cov_files = glob($cov_folder."/*.cov");
$cov_files[] = $cov;

$cov_files = array_unique(array_map("realpath", $cov_files));
if (count($cov_files)<$cov_min)
{
	generate_empty_cnv_file($out, $command, "", $ps_name, "Only ".count($cov_files)." coverage files found in folder '{$cov_folder}'", $tumor_only);
	trigger_error("CNV calling skipped. Only ".count($cov_files)." coverage files found in folder '$cov_folder'. At least {$cov_min} files are needed!", E_USER_ERROR);
}

//select coverage files of most similar samples
$corr_start = microtime(true);
$mean_correlation = 0.0;
{
	//sort coverage files
	sort($cov_files);
	
	//create target region without polymorphic regions
	$poly_merged = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedAdd", "-in {$repository_basedir}/data/misc/af_genomes_imgag.bed {$repository_basedir}/data/misc/centromer_telomer_hg19.bed -out {$poly_merged}", true);
	$roi_poly = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedIntersect", "-in {$bed} -in2 {$poly_merged} -out {$roi_poly} -mode in", true);
	$roi_nonpoly = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedSubtract", "-in {$bed} -in2 {$roi_poly} -out {$roi_nonpoly}", true);
	
	//determine which rows of coverage profiles to use
	$rows_to_use = determine_rows_to_use($cov, $roi_nonpoly);
	
	//determine correlation of each sample/cov file
	$file2corr = array();
	$cov1 = null;
	$cov2 = null;
	load_coverage_profile($cov, $rows_to_use, $cov1);
	foreach($cov_files as $cov_file)
	{
		load_coverage_profile($cov_file, $rows_to_use, $cov2);

		$corr = array();
		foreach($cov1 as $chr => $profile1)
		{
			$corr[] = number_format(correlation($profile1, $cov2[$chr]), 3);
		}
		$file2corr[$cov_file] = number_format(median($corr), 3);
	}
		
	//sort by correlation
	arsort($file2corr);
	
	//selects best n
	$file2corr = array_slice($file2corr, 0, $cov_max);
	$add_info = array();
	foreach ($file2corr as $f => $c)
	{
		$mean_correlation += $c;
		$add_info[] = "$f ($c)";
	}
	$mean_correlation = number_format($mean_correlation/count($file2corr), 3);
	
	$parser->log("Mean correlation to reference samples is {$mean_correlation}");
	$parser->log("Selected the following files as reference samples (correlation):", $add_info);
	
	$cov_files = array_keys($file2corr);
}
$parser->log("Execution time of determining reference samples: ".time_readable(microtime(true) - $corr_start));

//merge coverage files to one file
$tmp = $parser->tempFile(".txt");
sort($cov_files);
file_put_contents($tmp, implode("\n", $cov_files));
$cov_merged = $parser->tempFile(".cov");
$parser->exec(get_path("ngs-bits")."TsvMerge", "-in $tmp -cols chr,start,end -simple -out {$cov_merged}", true);

//collect off-target files for coverage files
$use_off_target = false;
if($tumor_only)
{
	if(isset($bed_off) && isset($cov_off))
	{
		if($cov_folder_off == "auto")
		{
			$cov_folder_off = "{$cov_folder}_off_target";
		}
		if(!is_dir($cov_folder_off))
		{
			trigger_error("CNV calling skipped. Off-target Coverage files folder cov_folder_n_off '$cov_folder_off' does not exist!", E_USER_ERROR);
		}

		//create file with off-target coverages
		$merged_cov_off = $parser->tempFile(".txt");
		$sample_ids = array();
		foreach($cov_files as $cov_name)
		{
			$sample_ids[] = basename($cov_name, ".cov");
		}
		$cov_paths_off = create_file_with_paths($cov_folder_off, realpath($cov_off), $sample_ids);
		$off_target_count = count($sample_ids);
		if ($off_target_count < $cov_off_min)
		{
			trigger_error("ClinCNV is called without off-target coverage files! A minimum of {$cov_off_min} files is necessary, {$off_target_count} given. Running ClinCNV for tumor-only analysis WITH off-target coverages is highly recommended!", E_USER_WARNING);
		}
		else
		{
			$parser->exec(get_path("ngs-bits")."/TsvMerge" , " -in $cov_paths_off -out {$merged_cov_off} -cols chr,start,end -simple",true, true);
			$use_off_target = true;
		}
	}
}

//execute ClinCNV (with workaround for hanging jobs)
$out_folder = $parser->tempFolder();
$args = [
"--normal {$cov_merged}",
"--bed {$bed}",
"--normalSample {$ps_name}",
"--out {$out_folder}",
"--numberOfThreads {$threads}",
"--par \"chrX:60001-2699520;chrX:154931044-155260560\"" //this is correct for hg19 only!
];

//analyzing a single tumor sample
if($tumor_only)
{
	$args[] = "--onlyTumor";
	$args[] = "--minimumPurity 30";
	$args[] = "--purityStep 5";
	$args[] = "--scoreS 150";
	$script_path = get_path('clincnv');
	$command_elements = explode(' ', $script_path);
	$script_path = end($command_elements);
	$args[] = "--folderWithScript {$script_path}";
	if($use_off_target)
	{
		$args[] = "--bedOfftarget $bed_off";
		$args[] = "--normalOfftarget $merged_cov_off";
	}
	if(is_dir($baf_folder)) $args[] = "--bafFolder {$baf_folder}";
}
else
{
	$args[] = "--maxNumGermCNVs {$max_cnvs}";
	$args[] = "--lengthG ".($regions-1); //lengthG actually gives the number of additional regions > subtract 1
	$args[] = "--scoreG 20";
}

//use_off_target is set for tumor_only with off target cov/bed files
if (!$skip_super_recall && !$use_off_target)
{
	$args[] = "--superRecall 3"; //superRecall will call down to one region and to log-likelihood 3
}

$parameters = implode(" ", $args);
$pid = getmypid();
$stdout_file = $parser->tempFile(".stdout", "megSAP_clincnv_pid{$pid}_");
$stderr_file = $parser->tempFile(".stderr", "megSAP_clincnv_pid{$pid}_");
$stdout = array();
$stderr = array();
$try_nr = 0;
$return = -1;
while($return!=0 && $try_nr < $max_tries)
{
	++$try_nr;
	
	//log
	$add_info = array();
	$add_info[] = "version    = ".$parser->extractVersion($command);
	$add_info[] = "parameters = $parameters";
	$parser->log("Calling external tool '$command' ({$try_nr}. try)", $add_info);

	$exec_start = microtime(true);
	$proc = proc_open($command." ".$parameters, array(1 => array('file',$stdout_file,'w'), 2 => array('file',$stderr_file,'w')), $pipes);
	while(true)
	{
		sleep(5);
		
		//check that clusters could be allocated
		$at_cluster_allocation = ends_with(trim(file_get_contents($stdout_file)), "\"START cluster allocation.\"");
		if ($at_cluster_allocation)
		{
			$sec_passed = time() - filemtime($stdout_file);
			if ($sec_passed>60)
			{
				$parser->log("Cluster allocation stuck in command '$command'. Stopping it and trying again...");
				proc_terminate($proc);
				$return = -999;
				break;
			}
		}
		
		$status  = proc_get_status($proc);
		if (!$status['running'])
		{
			$return = $status['exitcode'];
			break;
		}
	}
	
	//log output
	$stdout = explode("\n", rtrim(file_get_contents($stdout_file)));
	$stderr = explode("\n", rtrim(file_get_contents($stderr_file)));
	if (count($stdout)>0)
	{
		$parser->log("Stdout of '$command':", $stdout);
	}
	if (count($stderr)>0)
	{
		$parser->log("Stderr of '$command':", $stderr);
	}
}

//log error
$parser->log("Execution time of '$command': ".time_readable(microtime(true) - $exec_start));
if ($return!=0)
{
	$parser->toStderr($stdout);
	$parser->toStderr($stderr);
	trigger_error("Call of external tool '$command' returned error code '$return'.", E_USER_ERROR);
}

$clinCNV_result_folder = "normal";
if($tumor_only)
{
	$clinCNV_result_folder="tumorOnly";
}

if(file_exists("{$out_folder}/{$clinCNV_result_folder}/{$ps_name}/{$ps_name}_cnvs.tsv"))
{
	//sort and extract sample data from output folder
	$parser->exec(get_path("ngs-bits")."/BedSort","-in {$out_folder}/{$clinCNV_result_folder}/{$ps_name}/{$ps_name}_cnvs.tsv -out $out",true);
	$parser->copyFile("{$out_folder}/{$clinCNV_result_folder}/{$ps_name}/{$ps_name}_cov.seg", substr($out, 0, -4).".seg");
	$parser->copyFile("{$out_folder}/{$clinCNV_result_folder}/{$ps_name}/{$ps_name}_cnvs.seg", substr($out, 0, -4)."_cnvs.seg");

	//add high-quality CNV count to header
	$cnv_calls = file($out);
	$hq_cnvs = 0;
	$analysistype = 0;
	$i = 0;
	foreach($cnv_calls as $line)
	{
		if(starts_with($line, "##ANALYSISTYPE"))
		{
			$analysistype = $i;
		}
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		
		list($c, $s, $e, $tmp_cn, $ll) = explode("\t", $line);
		if ($ll>=20)
		{
			++$hq_cnvs;
		}	
		$i++;
	}
	array_splice($cnv_calls, 3, 0, array("##high-quality cnvs: {$hq_cnvs}\n", "##mean correlation to reference samples: {$mean_correlation}\n"));
	if($tumor_only)
	{
		array_splice($cnv_calls, $analysistype, 1, array("##ANALYSISTYPE=CLINCNV_TUMOR_ONLY\n"));
	}

	file_put_contents($out, $cnv_calls);
}
else
{
	//ClinVAR did not generate CNV file
	generate_empty_cnv_file($out, $command, $stdout, $ps_name,$stderr, $tumor_only);
}

?>

