<?php 
/** 
	@page vc_cnvhunter

	@todo check if filtering for germline CN-polymorphisms improves somatic output
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
$args[] = "-cnp_file ".repository_basedir()."/data/dbs/CNPs/copy_number_map_strict.bed";
$args[] = "-annotate ".repository_basedir()."/data/gene_lists/genes.bed ".repository_basedir()."/data/gene_lists/dosage_sensitive_disease_genes.bed {$data_folder}/dbs/OMIM/omim.bed";
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
		
		//check that at least half of the regions have a significantly higher number z-score in tumor than in normal
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
	
	//annotate marker positions, bed file needs to be sorted
	$tmp_file = $parser->tempFile("_sorted.bed");
	$parser->exec(get_path("ngs-bits")."BedSort", "-in ".$sys['target_file']." -out ".$tmp_file." -uniq", true);
	$target = Matrix::fromTSV($tmp_file);
	$new_col1 = array_fill(0,$cnvs_somatic->rows(),$target->rows());
	$new_col2 = array_fill(0,$cnvs_somatic->rows(),0);
	for($i=0;$i<$target->rows();++$i)
	{
		$row = $target->getRow($i);
		for($j=0;$j<$cnvs_somatic->rows();++$j)
		{
			$r = $cnvs_somatic->getRow($j);
			if($row[0]!=$r[0])	continue;
			if(!range_overlap($row[1],$row[2],$r[1],$r[2]))	continue;
			
			if($i<$new_col1[$j])	$new_col1[$j] = $i+1;
			if($i>$new_col2[$j])	$new_col2[$j] = $i+1;
		}
	}
	$cnvs_somatic->addCol($new_col1,"first_marker");
	$cnvs_somatic->addCol($new_col2,"last_marker");
	$idx_first = $cnvs_somatic->getColumnIndex("first_marker");
	$idx_last = $cnvs_somatic->getColumnIndex("last_marker");
	$idx_start = $cnvs_somatic->getColumnIndex("start");
	$idx_end = $cnvs_somatic->getColumnIndex("end");
	$tmp_matrix = new Matrix();
	$tmp_matrix->setComments($cnvs_somatic->getComments());
	$tmp_matrix->setHeaders($cnvs_somatic->getHeaders());
	
	//merge somatic regions
	if($cnvs_somatic->rows()>1)
	{
		$tmp_matrix->addRow($cnvs_somatic->getRow(0));
		for($i=1;$i<$cnvs_somatic->rows();++$i)
		{
			$current_region = $cnvs_somatic->getRow($i);
			$merged_region = $current_region;
			
			$merged_cn = "n/a";
			$tmp_median = median(explode(",",$merged_region[6]));
			if($tmp_median>2 && $tmp_median<4)
			{
				$merged_cn = "AMP";
			}
			else if($tmp_median>=4)
			{
				$merged_cn = "HIGH_AMP";
			}
			else if($tmp_median<2 && $tmp_median>0)	$merged_cn = "LOSS";
			else if($tmp_median<=0)	$merged_cn = "HIGH_LOSS";
			
			//merge regions only 10 % of markers / basepairs apart
			$max_dist = 0.1;
			if(contains("AMP",$merged_cn))	$max_dist = 0.15;

			//compare to previous regions
			$merge = true;
			$row_index = $tmp_matrix->rows()-1;
			while($merge && $row_index>=0)
			{
				$previous_region = $tmp_matrix->getRow($row_index);
				
				$previous_cn = "n/a";
				$tmp_median = median(explode(",",$previous_region[6]));
				if($tmp_median>2 && $tmp_median<4)	$previous_cn = "AMP";
				else if($tmp_median>=4)	$previous_cn = "HIGH_AMP";
				else if($tmp_median<2 && $tmp_median>0)	$previous_cn = "LOSS";
				else if($tmp_median<=0)	$previous_cn = "HIGH_LOSS";
				
				$allowed_markers = max(($merged_region[$idx_last]-$merged_region[$idx_first])*$max_dist,($previous_region[$idx_last]-$previous_region[$idx_first])*$max_dist);
				$allowed_dist = max(($merged_region[$idx_end]-$merged_region[$idx_start])*$max_dist,($previous_region[$idx_end]-$previous_region[$idx_start])*$max_dist);

				//check same type
				if($merged_region[0]!=$previous_region[0])	$merge = false;
				if($merged_cn!=$previous_cn)	$merge = false;
				if(($merged_region[$idx_first] - $previous_region[$idx_last]) > $allowed_markers && ($merged_region[$idx_start] - $previous_region[$idx_end]) >  $allowed_dist)	$merge = false;
				
				if($merge)
				{
					$merged_region[1] = $previous_region[1];
					$merged_region[4] += $previous_region[4];
					$merged_region[5] += $previous_region[5];
					$merged_region[6] = $previous_region[6].",".$merged_region[6];
					$merged_region[7] = $previous_region[7].",".$merged_region[7];
					$merged_region[8] = $previous_region[8].",".$merged_region[8];
					$merged_region[9] = $previous_region[9].",".$merged_region[9];
					$merged_region[10] = $previous_region[10].",".$merged_region[10];
					$merged_region[11] = $previous_region[11];
					
					// remove previous region
					$tmp_matrix->removeRow($row_index);
					--$row_index;
				}
			}
			
			$tmp_matrix->addRow($merged_region);
		}
	}
	$cnvs_somatic = $tmp_matrix;
	
	//fix region information
	//@todo use seg file to update marker count / zscores / regions
	$tmp_matrix = new Matrix();
	$tmp_matrix->setComments($cnvs_somatic->getComments());
	$tmp_matrix->setHeaders($cnvs_somatic->getHeaders());
	for($i=0;$i<$cnvs_somatic->rows();++$i)
	{
		$current_region = $cnvs_somatic->getRow($i);
	
		$count_regions = 0;
		$region_copy_numbers = array();
		$region_zscores = array();
		$region_coordinates = array();

		for($j=0;$j<$reduced_tumor->rows();++$j)
		{
			$rt = $reduced_tumor->getRow($j);

			list($chr,$region) = explode(":",$rt[1]);
			list($start,$end) = explode("-",$region);
			
			if($chr!=$current_region[0])	continue;
			if(!range_overlap($start,$end,$current_region[1],$current_region[2]))	continue;

			++$count_regions;
			$region_coordinates[] = $rt[1];
			$region_copy_numbers[] = $rt[2];
			$region_zscores[] = $rt[3];
		}
		
		$current_region[4] = $current_region[2]-$current_region[1]+1;
		$current_region[5] = $count_regions;
		$current_region[6] = implode(",",$region_copy_numbers);
		$current_region[7] = implode(",",$region_zscores);
		$current_region[8] = "no_data";	// information missing
		$current_region[9] = implode(",",$region_coordinates);
		$current_region[10] = "no_data";	// genes missing
		
		$tmp_matrix->addRow($current_region);
	}
	trigger_error("Columns GENES and AF may contain all regions due to the merging process.",E_USER_WARNING);
	$cnvs_somatic = $tmp_matrix;

	// anotate genes and replace old gene column
	$tmp_file1 = $temp_folder."/cnvs_tumor.tsv";
	$tmp_file2 = $temp_folder."/cnvs_tumor_ann.tsv";
	$idx_genes = $cnvs_somatic->getColumnIndex("genes");
	$cnvs_somatic->toTSV($tmp_file1);
	$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $tmp_file1 -out $tmp_file2 -extend 5", true);
	$tmp_somatic = Matrix::fromTSV($tmp_file2);
	$idx_last = $tmp_somatic->cols()-1;
	for($i=0;$i<$cnvs_somatic->rows();++$i)
	{
		if(array_slice($cnvs_somatic->getRow($i),0,3) != array_slice($tmp_somatic->getRow($i),0,3))	trigger_error("CNV files out of phase.",E_USER_ERROR);
		$cnvs_somatic->set($i,$idx_genes,str_replace(" ","",$tmp_somatic->get($i,$idx_last)));
	}
	$tmp_file3 = $temp_folder."/cnvs_tumor_ann_col.tsv";
	$cnvs_somatic->toTSV($tmp_file3);
	
	//decide if focal, chromosome_arm, chromosome and add column cnv_type
	//parse ucsc chromosome file
	$chr = "";
	$start = 0;
	$end = 0;
	$arm = "";
	$chrom_arms = array();
	$matrix = Matrix::fromTSV(repository_basedir()."/data/dbs/UCSC/cytoBand.txt");
	for($i=0;$i<$matrix->rows();++$i)
	{
		$row = $matrix->getRow($i);
		
		list($tmp_chr,$tmp_start,$tmp_end,$tmp_arm) = $row;
		++$tmp_start;	//0-based to 1-based
		$tmp_arm = $tmp_arm[0];
		
		if(empty($chr))	list($chr,$start,$end,$arm) = array($tmp_chr,$tmp_start,$tmp_end,$tmp_arm);
		else if($chr!=$tmp_chr || $arm!=$tmp_arm)
		{
			$chrom_arms[] = array($chr,$start,$end,$arm);
			list($chr,$start,$end,$arm) = array($tmp_chr,$tmp_start,$tmp_end,$tmp_arm);
		}
		else
		{
			$end = $tmp_end;
		}
	}
	if(!empty($chr))	$chrom_arms[] = array($chr,$start,$end,$arm);
	
	$new_col = array();
	$idx_genes = $cnvs_somatic->getColumnIndex("genes");
	for($i=0;$i<$cnvs_somatic->rows();++$i)
	{
		$row = $cnvs_somatic->getRow($i);
		
		$overlap = array();
		for($j=0;$j<count($chrom_arms);++$j)
		{
			$r = $chrom_arms[$j];
			if($row[0]!=$r[0])	continue;
			
			$intersect = range_intersect($row[1],$row[2],$r[1],$r[2]);
			$size = $r[2]-$r[1];
			
			// overlap - arm, fraction
			if(isset($overlap[$r[3]]))	trigger_error("This should not happen.",E_USER_ERROR);
			$overlap[$r[3]] = ($intersect[1]-$intersect[0])/$size;
		}
		
		if(count($overlap)>2 || count($overlap)<1)	trigger_error("This should not happen.",E_USER_ERROR);
		
		$p_arm = $overlap["p"];
		$q_arm = $overlap["q"];
		$chr = ($p_arm+$q_arm)/2;
		$acrocentric = false;
		if(in_array($row[0],array("chr13","chr14","chr15","chr21","chr22")))	$acrocentric = true;
		$tmp_type = "unknown";
		if($p_arm<=0.25 && $q_arm<=0.25)
		{
			$tmp_type = "focal";
			$count_genes = array_filter(explode(",",$row[$idx_genes]));
			if(count($count_genes)==0)	$tmp_type = "no gene";
			if(count($count_genes)>3)	$tmp_type = "cluster";
		}
		else if($p_arm>0.25 && $q_arm<=0.25)
		{
			$tmp_type = "partial p-arm";
			if($p_arm>0.75)	$tmp_type = "p-arm";
		}
		else if($q_arm>0.25 && $p_arm<=0.25)
		{
			$tmp_type = "partial q-arm";
			if($q_arm>0.75)	$tmp_type = "q-arm";
			if($q_arm>0.25 && $acrocentric)	$tmp_type = "partial chromosome";
			if($q_arm>0.75 && $acrocentric)	$tmp_type = "chromosome";
		}
		else if($p_arm>0.25 && $q_arm>0.25)
		{
			$tmp_type = "partial chromosome";
			if(($q_arm+$p_arm)/2>0.75) $tmp_type = "chromosome";
		}
				
		$new_col[] = $tmp_type;
	}
	$cnvs_somatic->insertCol(5,$new_col, "cnv_type","Type of CNV: focal (< 25 % of chrom. arm, < 3 genes), cluster (focal, > 3 genes), (partial) p/q-arm, (partial) chromosome.");

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
	$parser->moveFile($temp_folder."/cnvs.seg", substr($out, 0, -4).".seg");
}

?>