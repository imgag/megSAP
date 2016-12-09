<?php 
/** 
	@page vc_cnvhunter
	@todo somatic mode: check if overlapping regions are of same copy number state / also consider size of the overlap (in somatic mode - CS)
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_cnvhunter", "\$Rev: 911 $", "Wrapper for CnvHunter tool.");
$parser->addString("cov", "COV file of sample, e.g. 'GS123456_01.cov' (tumor COV file for somatic CNV calling).", false);
$parser->addOutfile("out", "Output CNV file in TSV format.", false);
//optional
$parser->addInfile("system", "Processing system INI file (determined from 'cov' sample name by default).", true);
$parser->addString("n_cov", "Normal COV file for somatic CNV calling.", true, null);
$parser->addString("cov_folder", "Folder with all coverage files (if different from [data_folder]/coverage/[system_short_name]/.", true, "auto");
$parser->addInt("min_reg", "Minimum number of subsequent regions for output CNV events.", true, 1);
$parser->addInt("n", "Number of (most similar) samples to consider.", true, 20);
$parser->addFloat("min_corr", "Minimum reference correlation for samples.", true, 0.8);
$parser->addOutfile("seg", "Output CNV file in SEG format (for IGV).", true);
extract($parser->parse($argv));

//init
$somatic = isset($n_cov) && !is_null($n_cov);
$psid1 = basename($cov,".cov");
$psid2 = $somatic ? basename($n_cov,".cov") : null;
$sys = load_system($system, $psid1);

//get coverage files (background)
if($cov_folder=="auto")
{
	$cov_folder = get_path("data_folder")."/coverage/".$sys['name_short'];
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
if ($somatic) $cov_files[] = $n_cov;

//convert paths to absolute paths and remove duplicates
$tmp_cov_files = array_unique(array_map("realpath", $cov_files));
foreach($tmp_cov_files as $key => $tcf)	if($tcf===FALSE)	trigger_error("Could not find coverage file ".$cov_files[$key],E_USER_ERROR);
$cov_files = $tmp_cov_files;

//run cnvhunter on all samples
$args = array();
$args[] = "-sam_min_corr $min_corr";
if($sys['type']=="WGS")
{
	$args[] = "-reg_min_cov 0.5";
	$args[] = "-sam_min_depth 0.5";
}
if(isset($seg)) $args[] = "-seg $seg";
$args[] = "-anno";
$args[] = "-n $n";
$temp = $parser->tempFile("_cnvhunter.tsv");
$base = substr($temp, 0, -4);
$parser->exec(get_path("ngs-bits")."CnvHunter", "-in ".implode(" ",$cov_files)." -out $temp ".implode(" ", $args), true);

//filter results for given processed sample(s)
$cnvs = Matrix::fromTSV($temp);
$cnvs_filtered = new Matrix();
$cnvs_filtered->setHeaders($cnvs->getHeaders());
for($i=0;$i<$cnvs->rows();++$i)
{
	$row = $cnvs->getRow($i);
	list($chr, $start, $end, $sample, $size, $num_reg) = $row;
	if ($num_reg<$min_reg) continue;
	if($sample==$psid1 || ($somatic && $sample==$psid2))
	{	
		$cnvs_filtered->addRow($row);
	}
}

//filter for differential CNVs
if($somatic)
{
	$cnvs_somatic = new Matrix();
	$cnvs_somatic->setHeaders($cnvs_filtered->getHeaders());
	for($i=0;$i<$cnvs_filtered->rows();++$i)
	{
		$row = $cnvs_filtered->getRow($i);
		list($chr, $start, $end, $sample, $size, $num_reg) = $row;
		
		if($sample!=$psid1) continue; //skip normal sample entries
		
		$overlap = false;
		for($j=0;$j<$cnvs_filtered->rows();++$j)
		{
			$r = $cnvs_filtered->getRow($j);
			if($r[3]!=$psid2) continue;
			if ($r[0]!=$chr) continue; //same chromosome
			if (range_overlap($start, $end, $r[1], $r[2]))
			{
				$overlap = true;
				break;
			}
		}

		//filter somatic variants
		if(!$overlap)
		{
			$cnvs_somatic->addRow($row);
		}
	}
	$cnvs_filtered = $cnvs_somatic;
}

//extract sample info and statistics of all samples
$qc_problems = false;
$hits = array();
$median_correl = array();
$median_cnvs = array();
$sample_info = file($base."_samples.tsv");
foreach($sample_info as $line)
{
	list($sample, , , $ref_correl, , $cnvs, $qc_info) = explode("\t", nl_trim($line));
	if ($sample==$psid1 || ($somatic && $sample==$psid2))
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
	copy2($base.".seg", substr($out, 0, -4).".seg");
}

?>