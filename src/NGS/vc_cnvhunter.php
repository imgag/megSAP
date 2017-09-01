<?php 
/** 
	@page vc_cnvhunter
	
	@todo check if filtering for germline CN-polymorphisms improved somatic output
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_cnvhunter", "Wrapper for CnvHunter tool.");
$parser->addString("cov", "COV file of sample, e.g. 'GS123456_01.cov' (tumor COV file for somatic CNV calling).", false);
$parser->addOutfile("out", "Output CNV file in TSV format.", false);
//optional
$parser->addInfile("system", "Processing system INI file (determined from 'cov' sample name by default).", true);
$parser->addString("n_cov", "Normal COV file for somatic CNV calling.", true, null);
$parser->addString("debug", "Folder for debug information.", true, null);
$parser->addString("cov_folder", "Folder with all coverage files (if different from [data_folder]/coverage/[system_short_name]/.", true, "auto");
$parser->addInt("n", "Number of (most similar) samples to consider.", true, 30);
$parser->addFloat("min_corr", "Minimum reference correlation for samples.", true, 0.8);
$parser->addFloat("min_z", "Minimum z-score to call a CNV.", true, 4.0);
$parser->addString("seg", "Sample name for which a SEG file is generated.", true);
extract($parser->parse($argv));

//init
$somatic = isset($n_cov) && !is_null($n_cov);
$psid1 = basename($cov,".cov");
$psid2 = $somatic ? basename($n_cov,".cov") : null;
$sys = load_system($system, $psid1);
$data_folder = get_path("data_folder");

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
if ($somatic) $cov_files[] = $n_cov;

//convert paths to absolute paths and remove duplicates
$tmp_cov_files = array_unique(array_map("realpath", $cov_files));
foreach($tmp_cov_files as $key => $tcf)	if($tcf===FALSE)	trigger_error("Could not find coverage file ".$cov_files[$key],E_USER_ERROR);
$cov_files = $tmp_cov_files;

//run cnvhunter on all samples
$args = array();
$args[] = "-min_z $min_z";
$args[] = "-sam_min_corr $min_corr";
if($sys['type']=="WGS")
{
	$args[] = "-reg_min_cov 0.5";
	$args[] = "-sam_min_depth 0.5";
}
if(isset($seg) && is_null($n_cov)) $args[] = "-seg $seg";
$args[] = "-n $n";
if ($somatic) $args[] = "-debug ALL";
$temp_folder = !empty($debug) ? $debug : $parser->tempFolder();
if(!is_dir($temp_folder))	mkdir ($temp_folder);
$args[] = "-cnp_file {$data_folder}/dbs/CNPs/copy_number_map_strict.bed";
$args[] = "-annotate {$data_folder}/gene_lists/genes.bed {$data_folder}/gene_lists/dosage_sensitive_disease_genes.bed {$data_folder}/dbs/OMIM/omim.bed";
$parser->exec(get_path("ngs-bits")."CnvHunter", "-in ".implode(" ",$cov_files)." -out ".$temp_folder."/cnvs.tsv ".implode(" ", $args), true);

// filter results for given processed sample(s)
$cnvs_unfiltered = Matrix::fromTSV($temp_folder."/cnvs.tsv");
$cnvs_filtered = new Matrix();
$cnvs_filtered->setHeaders($cnvs_unfiltered->getHeaders());
for($i=0; $i<$cnvs_unfiltered->rows(); ++$i)
{
	$line = $cnvs_unfiltered->getRow($i);
	list($chr, $start, $end, $sample, $size, $num_reg) = $line;

	if($sample==$psid1)
	{
		$cnvs_filtered->addRow($line);
	}
}

if($somatic)
{	
	// filter for tumor/normal debug data
	$reduced_tumor = new Matrix();
	$reduced_normal = new Matrix();
	$z_normal = array();
	$handle = fopen($temp_folder."/cnvs_debug.tsv", "r", "r");
	while(!feof($handle))
	{
		$line = nl_trim(fgets($handle));
		if(empty($line))	continue;
		if($line[0]=="#")	continue;
		
		$row = explode("\t",$line);
		if($row[0]==$psid1)
		{
			$reduced_tumor->addRow($row);
		}
		if($row[0]==$psid2)
		{
			$z_normal[$row[1]] = $row[3];
			$reduced_normal->addRow($row);
		}
	}

	//generate differential SEG file
	if(isset($seg))
	{
		$seg_file = array();
		$seg_file[] = "#type=GENE_EXPRESSION";
		$seg_file[] = "#track graphtype=points name=\"".$psid1."-".$psid2." CN z-score\" midRange=-4:4 color=0,0,255 altColor=255,0,0 viewLimits=-10:10 maxHeightPixels=80:80:80";
		$seg_file[] = "ID\tchr\tstart\tend\ttum z-score\tnor z-score\ttum ndoc\tnor ndoc\ttum ref_ndoc\tnor ref_ndoc\tcopy-number\tz-score";

		//write valid region details
		for($i=0;$i<$reduced_tumor->rows();++$i)
		{
			$row_t = $reduced_tumor->getRow($i);
			$row_n = $reduced_normal->getRow($i);
			if($row_t[1]!=$row_n[1])	trigger_error("Could not generate somatic SEG file. Regions in input coverage files are unsorted.",E_USER_ERROR);
			list($chr,$region) = explode(":",$row_t[1]);
			list($start,$end) = explode("-",$region);
			
			$seg_file[] = $psid1."-".$psid2."\t".$chr."\t".$start."\t".$end."\t".$row_t[3]."\t".$row_n[3]."\t".$row_t[4]."\t".$row_n[4]."\t".$row_t[5]."\t".$row_n[5]."\t".$row_t[2]."\t".$row_t[3];
		}
		
		file_put_contents($temp_folder."/cnvs.seg",implode("\n",$seg_file));
	}

	//filter somatic variants
	$cnvs_somatic = new Matrix();
	$headers = array_slice($cnvs_filtered->getHeaders(), 0, 10);
	$headers[] = "genes";
	$cnvs_somatic->setHeaders($headers);
	for($i=0;$i<$cnvs_filtered->rows();++$i)
	{
		$row = $cnvs_filtered->getRow($i);
		list($chr, $start, $end, $sample, $size, $num_reg, $reg_cns, $reg_zs, $reg_afs, $reg_coords, $overlap_cnp_region, $genes) = $row;
		
		$reg_cns = explode(",", $reg_cns);
		$reg_zs = explode(",", $reg_zs);
		$reg_coords = explode(",", $reg_coords);
		$reg_afs = explode(",", $reg_afs);

		//check that z-score of tumor and normal sample differ by at least factor 2
		$z_diff_pass = array();
		for ($r=0; $r<count($reg_coords); ++$r)
		{
			$reg = $reg_coords[$r];
			$z_t = $reg_zs[$r];
			$z_n = $z_normal[$reg];
			
			$pass = true;
			if ($z_t>0 && $z_n>0 && $z_t<2*$z_n)
			{
				$pass = false;
			}
			if ($z_t<0 && $z_n<0 && $z_t>2*$z_n)
			{
				$pass = false;
			}
			$z_diff_pass[] = $pass;
			//print "$reg\t$z_t\t$z_n\t$pass\n";
		}
		
		//trim front
		$start_index = null;
		for ($r=0; $r<count($reg_coords); ++$r)
		{
			$z_t = $reg_zs[$r];
			if ($z_diff_pass && ($z_t>3.0 || $z_t<-3.0))
			{
				$start_index = $r;
				break;
			}
		}
		
		//trim end
		$end_index = null;
		for ($r=count($reg_coords)-1; $r>=0; --$r)
		{
			$z_t = $reg_zs[$r];
			if ($z_diff_pass && ($z_t>3.0 || $z_t<-3.0))
			{
				$end_index = $r;
				break;
			}
		}
		
		//check number of regions (at least 3)
		$reg_count = $end_index - $start_index + 1;
		if (is_null($start_index) || is_null($end_index) || $reg_count < 3)
		{
			//print "FAIL: num_reg=$reg_count\n";
			continue;
		}
		
		//check that at least half of the regions have a significantly hiher number z-score in tumor than in normal
		$z_diff_pass_count = array_sum(array_slice($z_diff_pass, $start_index, $reg_count));
		if ($z_diff_pass_count < 0.5 * $reg_count)
		{		
			//print "FAIL: z_diff_pass_count=$z_diff_pass_count (of $reg_count)\n";
			continue;
		}
		
		//line passed => construct new line based on new boundaries
		$tmp  = explode(":", strtr($reg_coords[$start_index], "-", ":"));
		$start = $tmp[1];
		$tmp  = explode(":", strtr($reg_coords[$end_index], "-", ":"));
		$end = $tmp[2];
		$reg_cns = implode(",", array_slice($reg_cns, $start_index, $reg_count));
		$reg_zs = implode(",", array_slice($reg_zs, $start_index, $reg_count));
		$reg_coords = implode(",", array_slice($reg_coords, $start_index, $reg_count));
		$reg_afs = implode(",", array_slice($reg_afs, $start_index, $reg_count));
		$cnvs_somatic->addRow(array($chr, $start, $end, $sample, $end-$start+1, $reg_count, $reg_cns, $reg_zs, $reg_afs, $reg_coords, $genes));
	}
	$cnvs_filtered = $cnvs_somatic;
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
	copy2($temp_folder."/cnvs.seg", substr($out, 0, -4).".seg");
}

?>