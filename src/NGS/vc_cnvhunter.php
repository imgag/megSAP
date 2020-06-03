<?php 
/** 
	@page vc_cnvhunter
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_cnvhunter", "Wrapper for CnvHunter tool.");
$parser->addString("cov", "COV file of sample, e.g. 'GS123456_01.cov'.", false);
$parser->addOutfile("out", "Output CNV file in TSV format.", false);
//optional
$parser->addInfile("system", "Processing system INI file (determined from 'cov' sample name by default).", true);
$parser->addString("debug", "Folder for debug information.", true, null);
$parser->addString("cov_folder", "Folder with all coverage files (if different from [data_folder]/coverage/[system_short_name]/.", true, "auto");
$parser->addInt("n", "Number of (most similar) samples to consider.", true, 30);
$parser->addFloat("min_corr", "Minimum reference correlation for samples.", true, 0.8);
$parser->addFloat("min_z", "Minimum z-score to call a CNV.", true, 4.0);
$parser->addString("seg", "Sample name for which a SEG file is generated.", true);
extract($parser->parse($argv));

//init
$psid = basename($cov,".cov");
$sys = load_system($system, $psid);
$data_folder = get_path("data_folder");
$repository_basedir = repository_basedir();
$ngsbits = get_path("ngs-bits");

//get coverage files (background)
if($cov_folder=="auto")
{
	$cov_folder = "{$data_folder}/coverage/".$sys['name_short'];
}
if(!is_dir($cov_folder))
{
	trigger_error("CNV calling skipped. Coverage files folder '$cov_folder' does not exist!", E_USER_WARNING);
	exit(0);
}
$cov_files = glob($cov_folder."/*.cov");
if (count($cov_files)<$n+1)
{
	trigger_error("CNV calling skipped. Only ".count($cov_files)." coverage files found in folder '$cov_folder'. At least ".($n+1)." files are needed!", E_USER_WARNING);
	exit(0);
}
$cov_files[] = $cov;

//convert paths to absolute paths, remove duplicates and store input file list
$tmp_cov_files = array_unique(array_map("realpath", $cov_files));
foreach($tmp_cov_files as $key => $tcf)	if($tcf===FALSE)	trigger_error("Could not find coverage file ".$cov_files[$key],E_USER_ERROR);
$cov_files = $tmp_cov_files;
$cov_list = $parser->tempFile("_cov_files.txt");
file_put_contents($cov_list, implode("\n", $cov_files));

//run cnvhunter on all samples
$args = array();
$args[] = "-in {$cov_list}";
$args[] = "-min_z $min_z";
$args[] = "-sam_min_corr $min_corr";
if($sys['type']=="WGS" || $sys['type']=="WGS (shallow)")
{
	$args[] = "-reg_min_cov 0.5";
	$args[] = "-sam_min_depth 0.5";
}
else
{
	$args[] = "-sam_min_depth 10";
}
if(isset($seg)) $args[] = "-seg $seg";
$args[] = "-n $n";
$temp_folder = !empty($debug) ? $debug : $parser->tempFolder();
if(!is_dir($temp_folder) || !is_writable($temp_folder))
{
	trigger_error("Temp folder '$temp_folder' not writable.", E_USER_ERROR);
}
$out_tmp = "{$temp_folder}/cnvs.tsv";
$args[] = "-out {$out_tmp}";
$args[] = "-cnp_file {$repository_basedir}/data/misc/af_genomes_imgag.bed";
$parser->exec("{$ngsbits}CnvHunter", implode(" ", $args), true);

// filter results for given processed sample(s)
$cnvs_unfiltered = Matrix::fromTSV($out_tmp);
$cnvs_filtered = new Matrix();
$cnvs_filtered->setHeaders($cnvs_unfiltered->getHeaders());
for($i=0; $i<$cnvs_unfiltered->rows(); ++$i)
{
	$line = $cnvs_unfiltered->getRow($i);
	list($chr, $start, $end, $sample, $size, $num_reg) = $line;

	if($sample==$psid)
	{
		$cnvs_filtered->addRow($line);
	}
}

//extract sample info and statistics of all samples
$qc_problems = false;
$hits = array();
$median_correl = array();
$median_cnvs = array();
$sample_info = file($temp_folder."/cnvs_samples.tsv");
foreach($sample_info as $line)
{
	list($sample, , $ref_correl, , $cnvs, $qc_info) = explode("\t", nl_trim($line));
	if ($sample==$psid)
	{
		$hits[] = array($sample, $ref_correl, $cnvs, $qc_info);
		if ($qc_info!="")
		{
			trigger_error("CNV calling of sample $sample not possible: $qc_info!", E_USER_WARNING);
			$qc_problems = true;
		}
	}
	else
	{
		$median_correl[] = $ref_correl;
		$median_cnvs[] = $cnvs;
	}
}

//add analysis infos to TSV header
$cnvs_filtered->addComment("#ANALYSISTYPE=CNVHUNTER_GERMLINE_SINGLE");
list($cnvhunter_version) = exec2("{$ngsbits}CnvHunter --version");
$cnvhunter_version = trim(substr($cnvhunter_version[0], 9));
$cnvs_filtered->addComment("#CnvHunter version: $cnvhunter_version");

//add sample info to TSV header
$median_correl = median($median_correl);
$median_cnvs = median($median_cnvs);
foreach($hits as $values)
{
	list($sample, $ref_correl, $cnvs, $qc_info) = $values;
	$cnvs_filtered->addComment("#$sample ref_correl: $ref_correl (median: $median_correl)");
	$cnvs_filtered->addComment("#$sample cnvs: $cnvs (median: $median_cnvs)");
	$cnvs_filtered->addComment("#$sample qc_errors: $qc_info");
}

//store output tsv file
$cnvs_filtered->toTSV($out);
if(isset($seg) && !$qc_problems)
{
	$parser->moveFile($temp_folder."/cnvs.seg", substr($out, 0, -4).".seg");
}

?>