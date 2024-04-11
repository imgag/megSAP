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
$parser->addInt("cov_min", "Minimum number of reference coverage files required for CNV analysis.", true, 10);
$parser->addInt("cov_max", "Maximum number of reference coverage files used for CNV analysis. This parameter is needed to keep run-time and RAM requirement manageable.", true, 150);
$parser->addInt("cov_compare_max", "Maximum number of coverage files to compare during similarity calculation. Only possible with NGSD support enabled.", true, 600);
$parser->addInt("max_cnvs", "Number of expected CNVs (~200 for WES and ~2000 for WGS).", true, 2000);
$parser->addInt("max_tries", "Maximum number of tries for calling ClinCNV (R parallelization sometimes breaks with no reason", true, 10);
$parser->addInt("regions", "Number of subsequent regions that must show a signal for a call.", true, 2);
$parser->addFlag("skip_super_recall", "Skip super-recall (down to one region and log-likelihood 3).");
$parser->addFlag("mosaic", "Detect additionally large mosaic regions (for WES, WGS and shallow WGS");
$parser->addEnum("gender", "Force gender used by ClinCNV", true, ['male', 'female', 'n/a'], 'n/a');

//tumor only flags (use for somatic with no normal sample)
$parser->addFlag("tumor_only", "Analyze tumor sample without a paired normal sample.");
$parser->addInfile("bed_off","Off-target bed file.",true); //s_dna
$parser->addInfile("cov_off","Off-target coverage file for normal sample",true);
$parser->addString("cov_folder_off", "Folder with off-target normal coverage files (if different from [data_folder]/coverage/[system_short_name]_off_target/.", true, "auto");
$parser->addString("baf_folder","Folder containing files with B-Allele frequencies.",true);
$parser->addString("cov_off_min","Minimum number of off-target coverage files required for CNV analysis.",true, 10);

extract($parser->parse($argv));

class cnv {
	public $chr;
	public $start;
	public $end;
} 

function cnv_is_in_list($cnv, $list)
{
	$ori_length = $cnv->end-$cnv->start;
	$substr_length = $ori_length;
	foreach($list as $region)
	{
		if($region->chr == $cnv->chr)
		{
			//left overlap
			if($cnv->end > $region->start && $cnv->end <= $region->end && $cnv->start < $region->start)
			{
				$substr_length -= $cnv->end-$region->start;
			}
			//right overlap
			else if($cnv->end > $region->end && $cnv->start < $region->end && $cnv->start >= $region->start)
			{
				$substr_length -= $region->end-$cnv->start;
			}
			//inside
			else if($cnv->end <= $region->end && $cnv->start >= $region->start)
			{
				return 100;
			}
			//outside
			else if($cnv->start < $region->start && $cnv->end > $region->end)
			{
				$substr_length -= $region->end-$region->start;
			}
		}
	}
	//if at least 50% of the regions map, don't report it as mosaic CNV
	if($substr_length <= 0.5*$ori_length) return true;
	return false;	
}

//filter mosaic CNVs by already found CNVs
function filterMosaicVariants($mosaic_out, $filter_regions)
{

	$mosaicism = false;

	//get possible mosaic cnvs
	$mosaic_h = fopen2($mosaic_out, "r");
	$new_lines = array();
	$length_kb_idx = 0;
    $cn_idx = 0;
	while(!feof($mosaic_h))
	{
		$line = fgets($mosaic_h);

		if(starts_with(trim($line), "#chr"))
        {
            $count = 0;
            $entries = explode("\t", $line);
            foreach($entries as $col)
            {
                if(trim($col) == "length_KB")
                {
                    $length_kb_idx = $count;
                }
                else if(trim($col) == "CN_change")
                {
                    $cn_idx = $count;
                }
                
                $count++;
			}
			
			$new_lines[] = $line;
        }
		else if(starts_with(trim($line), "#"))
		{
			$new_lines[] = $line;
		}
		else
		{
			$entries = explode("\t", $line);
			if(count($entries) >= max(3, $cn_idx, $length_kb_idx))
			{
				$m_cnv = new cnv();

				$m_cnv->chr = $entries[0]; 
				$m_cnv->start = $entries[1]; 
				$m_cnv->end = $entries[2]; 

				$mosaic_cn = $entries[$cn_idx];
				$mosaic_cn = doubleval($mosaic_cn);
				$mosaic_length_bases =  str_replace('.', '', $entries[$length_kb_idx]);
				$mosaic_length_bases = intval($mosaic_length_bases);

				if(!cnv_is_in_list($m_cnv, $filter_regions) && $mosaic_length_bases >= 500000 && 1<$mosaic_cn && $mosaic_cn<3)
				{
					$new_lines[] = $line;
					$mosaicism = true;
				}
			}
		}
	}
	fclose($mosaic_h);

	$mosaic_h = fopen2($mosaic_out, "w");
	foreach($new_lines as $line)
	{
		fwrite($mosaic_h, $line);
	}
	fclose($mosaic_h);

	return $mosaicism;
}

function detect_mosaicism()
{
	global $parser;
	global $out;
	global $repository_basedir;

	//run ClinCNV to detect large mosaic CNVs
	$tmp_folder = $parser->tempFolder();
	$mosaic_out = $tmp_folder."/mosaic.tsv";
	
	$mosaicism = TRUE;
	run_clincnv($mosaic_out, $mosaicism);

	//merge bed regions of polymorphism and CNVs
	$polymorphic="{$repository_basedir}/data/misc/af_genomes_imgag.bed";
	$combined_bed = $parser->tempFile(".bed");
	$filter_regions_bed = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedAdd", "-in {$polymorphic} {$out} -out {$combined_bed}", true);
	$parser->exec(get_path("ngs-bits")."BedMerge", "-in {$combined_bed} -out {$filter_regions_bed}", true);

	//store all those regions
    $regions_filter = array();
    $filter_h = fopen2($filter_regions_bed, "r");
    while(!feof($filter_h))
    {
        $line = trim(fgets($filter_h));

        if(starts_with($line, "#"))
        {
            continue;
        }
        else
        {
            $entries = explode("\t", $line);
            if(count($entries) >= 3)
            {
                
                $tmp_cnv = new cnv();
                $tmp_cnv->chr = $entries[0];
                $tmp_cnv->start = $entries[1];
                $tmp_cnv->end = $entries[2];

                $old_cnv = end($regions_filter);
                if(!empty($old_cnv) && ($tmp_cnv->start <= $old_cnv->end) && $tmp_cnv->chr == $old_cnv->chr)
                {
                    array_pop($regions_filter);
                    $old_cnv->end = $entries[2];
                    $regions_filter[] = $old_cnv;
                }
                else
                {				
                    $regions_filter[] = $tmp_cnv;
                }	
            }
        }
	}
	fclose($filter_h);

	filterMosaicVariants($mosaic_out, $regions_filter);

	//copy results to output folder
	$sample_cnv_name = substr($out,0,-4);
	$mosaic_file = $sample_cnv_name."_mosaic.tsv";
	if (file_exists($mosaic_out)) $parser->moveFile($mosaic_out, $mosaic_file);
}

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
function generate_empty_cnv_file($out, $command, $stdout, $ps_name, $error_messages, $tumor_only)
{
	//ClinCNV did not generate CNV file
	//generate file with basic header lines
	$cnv_output = fopen2($out, "w");
	if($tumor_only)
	{
		fwrite($cnv_output, "##ANALYSISTYPE=CLINCNV_TUMOR_ONLY\n");
	}
	else
	{
		fwrite($cnv_output, "##ANALYSISTYPE=CLINCNV_GERMLINE_SINGLE\n");
	}
	fwrite($cnv_output, "##GENOME_BUILD=GRCh38\n");
	if(preg_match('/^.*ClinCNV-([\d\.]*).*$/', $command, $matches))
	{
		fwrite($cnv_output, "##ClinCNV version: v{$matches[1]}\n");
	}
	if (!is_array($error_messages))
	{
		$error_messages = array($error_messages);
	}
	foreach($error_messages as $error_message)
	{
		fwrite($cnv_output, "##CNV calling skipped: {$error_message}.\n");
	}
	//search in CLinCNV stdout for lines indicating failed QC
	foreach($stdout as $line)
	{
		//if QC for input sample failed, add it to header line
		if(contains($line, "{$ps_name} did not pass QC"))
		{
			if(preg_match('/^.*\"(.*)\".*$/', $line, $matches))
			{
				fwrite($cnv_output, "##{$matches[1]}\n");
			}
		}
	}
	fwrite($cnv_output, "#chr\tstart\tend\tCN_change\tloglikelihood\tno_of_regions\tlength_KB\tpotential_AF\tgenes\n");

	fclose($cnv_output);
}

function run_clincnv($out, $mosaic=FALSE)
{
	global $tumor_only;
	global $use_off_target;
	global $bed_off;
	global $merged_cov_off;
	global $baf_folder;
	global $max_cnvs;
	global $regions;
	global $max_tries;
	global $command;
	global $parser;
	global $mean_correlation;
	global $ps_name;
	global $skip_super_recall;
	global $cov_merged;
	global $bed;
	global $threads;
	global $gender;

	$out_folder = $parser->tempFolder();

	$args = [
		"--normal {$cov_merged}",
		"--bed {$bed}",
		"--normalSample {$ps_name}",
		"--out {$out_folder}",
		"--numberOfThreads {$threads}",
		"--par \"chrX:10001-2781479;chrX:155701383-156030895\"", //this is correct for hg38 only!
		"--hg38",
		"--noPlot",
		"--folderWithScript ".dirname(get_path("clincnv"))
		];
		
	if ($gender=="male") $args[] = "--sex M";
	if ($gender=="female") $args[] = "--sex F";

	//analyzing a single tumor sample
	if($tumor_only)
	{
		$args[] = "--onlyTumor";
		$args[] = "--minimumPurity 30";
		$args[] = "--purityStep 5";
		$args[] = "--scoreS 150";
		if($use_off_target)
		{
			$args[] = "--bedOfftarget $bed_off";
			$args[] = "--normalOfftarget $merged_cov_off";
		}
		if(is_dir($baf_folder)) $args[] = "--bafFolder {$baf_folder}";
	}
	else if($mosaic)
	{
		$args[] = "--maxNumGermCNVs 10";
		$args[] = "--lengthG 2"; //lengthG actually gives the number of additional regions > subtract 1
		$args[] = "--scoreG 500";
		$args[] = "--mosaicism";
	}
	else
	{
		$args[] = "--maxNumGermCNVs {$max_cnvs}";
		$args[] = "--lengthG ".($regions-1); //lengthG actually gives the number of additional regions > subtract 1
		$args[] = "--scoreG 20";
		
	}

	//use_off_target is set for tumor_only with off target cov/bed files
	if (!$skip_super_recall && !$use_off_target && !$mosaic)
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

	if(!file_exists("{$out_folder}/{$clinCNV_result_folder}/{$ps_name}/{$ps_name}_cnvs.tsv"))
	{
		//ClinCNV did not generate CNV file
		generate_empty_cnv_file($out, $command, $stdout, $ps_name, $stderr, $tumor_only);
		return false;
	}
	
	//sort and extract sample data from output folder
	$parser->exec(get_path("ngs-bits")."/BedSort","-in {$out_folder}/{$clinCNV_result_folder}/{$ps_name}/{$ps_name}_cnvs.tsv -out $out",true);
	$parser->copyFile("{$out_folder}/{$clinCNV_result_folder}/{$ps_name}/{$ps_name}_cov.seg", substr($out, 0, -4).".seg");
	$parser->copyFile("{$out_folder}/{$clinCNV_result_folder}/{$ps_name}/{$ps_name}_cnvs.seg", substr($out, 0, -4)."_cnvs.seg");
	
	//post-processing of calls:
	//- add analysis type, high-quality CNV count and mean correlation to reference samples to header
	//- remove call with log-likelihood below 0 (they are artefacts: ClinCNV first does the calling and the a correction of log-likelihoods. In very noisy regions that log-likelihood can be below zero after the correction)
	$cnv_calls = [];
	$hq_cnvs = 0;
	$analysistype = 0;
	$i = 0;
	$h = fopen2($out, "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="") continue;
		
		if ($line[0]=="#") //header line
		{
			//replace type header for somatic
			if($tumor_only && starts_with($line, "##ANALYSISTYPE"))
			{
				$line = "##ANALYSISTYPE=CLINCNV_TUMOR_ONLY";
			}
		}
		else //content line
		{
			list($c, $s, $e, $cn, $ll) = explode("\t", $line);
			
			if ($ll>=20) ++$hq_cnvs;
			
			if ($ll<0) continue;
		}
		
		$cnv_calls[] = $line."\n";
	}
	fclose($h);
	
	$add_header_lines = [];
	$add_header_lines[] = "##GENOME_BUILD=GRCh38\n";
	$add_header_lines[] = "##high-quality cnvs: {$hq_cnvs}\n";
	$add_header_lines[] = "##mean correlation to reference samples: {$mean_correlation}\n";
	array_splice($cnv_calls, 3, 0, $add_header_lines);

	file_put_contents($out, $cnv_calls);
	
	return true;
}


//init
$repository_basedir = repository_basedir();
$ps_name = basename($cov,".cov");
$command = get_path("rscript")." --vanilla ".get_path("clincnv");

//determine coverage files
$cov_files = glob($cov_folder."/*.cov");
$cov_files[] = $cov;
$cov_files = array_unique(array_map("realpath", $cov_files));
if (count($cov_files)<$cov_min)
{
	generate_empty_cnv_file($out, $command, [], $ps_name, "Only ".count($cov_files)." coverage files found in folder '{$cov_folder}'", $tumor_only);
	trigger_error("CNV calling skipped. Only ".count($cov_files)." coverage files found in folder '$cov_folder'. At least {$cov_min} files are needed!", E_USER_ERROR);
}

//if too many, select the samples processed roughly at the same time based on sequencing run date
$restrict_cov_file_number = (count($cov_files)> $cov_compare_max) && db_is_enabled("NGSD");
if ($restrict_cov_file_number)
{
	$parser->log("Restricting number of coverage files to $cov_compare_max (based on sequencing run date from NGSD)...");
	
	$db = DB::getInstance("NGSD", false);
	$run_name = $db->getValue("SELECT r.name FROM sequencing_run r, processed_sample ps WHERE ps.sequencing_run_id=r.id AND ps.id=".get_processed_sample_id($db, $ps_name, false), "");
	if ($run_name!="#00000")
	{
		$ps2date = [];
		$res = $db->executeQuery("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) as name, r.start_date as date FROM processed_sample as ps, sample as s, sequencing_run r WHERE ps.sample_id = s.id AND ps.sequencing_run_id=r.id");
		foreach($res as $row)
		{
			$ps = trim($row['name']);
			$date = trim($row['date']);
			if ($date=="") continue;
			
			$ps2date[$ps] = $date;
		}
		
		if (isset($ps2date[$ps_name]))
		{
			$ps_time = strtotime($ps2date[$ps_name]);
			$cov2timediff = [];
			foreach($cov_files as $cov_file)
			{
				$ps_cov = basename($cov_file, ".cov");
				$ps_date = isset($ps2date[$ps_cov]) ? $ps2date[$ps_cov] : "2000-01-01";
				$cov2timediff[$cov_file] = abs($ps_time - strtotime($ps_date));
			}
			asort($cov2timediff);
			$cov_files = array_keys($cov2timediff);
		}
		else
		{
			$parser->log("Notice: Could not restrict the number of coverage files: NGSD does not contain sample '$ps_name'!");
			$restrict_cov_file_number = false;
		}
	}
	else
	{
			$parser->log("Notice: Could not restrict the number of coverage files: Sequencing run is dummy run '#00000'!");
			$restrict_cov_file_number = false;
	}
}

//select coverage files of most similar samples
$corr_start = microtime(true);
$mean_correlation = 0.0;
{	
	//create target region without polymorphic regions
	$poly_merged = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedAdd", "-in {$repository_basedir}/data/misc/af_genomes_imgag.bed {$repository_basedir}/data/misc/centromer_telomer.bed -out {$poly_merged}", true);
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
		
		if ($restrict_cov_file_number && count($file2corr)>=$cov_compare_max) break;
	}
	$parser->log("Compared number of coverage files: ".count($file2corr));

	//selects best n by correlation
	arsort($file2corr);
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

//call CNVs
$cnvs_called=run_clincnv($out);

//call mosaic CNVs
if($mosaic && $cnvs_called && !$tumor_only)
{
	detect_mosaicism();
}

?>
