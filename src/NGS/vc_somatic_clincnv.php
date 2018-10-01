<?php 
/** 
	@page vc_somatic_clincnv

*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_somatic_clincnv", "Wrapper for CnvHunter tool.");
$parser->addString("t_id", "Processed sample id of tumor, e.g. 'GS123456_01'.", false);
$parser->addInfile("t_cov","Coverage file for tumor sample",false);
$parser->addInfile("n_cov","Coverage file for normal sample",false);
$parser->addOutFile("out", "Output file.",false);
//optional
$parser->addInfile("bed","Bed file for target region.",false);
$parser->addString("cohort_folder", "Folder that contains ClinCNV cohort data (same processing system).", true,"auto");
$parser->addString("n_id","Processed sample id of normal sample, e.g. 'GS123456_01'. Using normal sample ID from NGSD by default.",true,"auto");
$parser->addString("cov_folder_n", "Folder with all coverage files (if different from [data_folder]/coverage/[system_short_name]/.", true, "auto");
$parser->addString("cov_folder_t", "Folder with all tumor coverage files (if different from [data_folder]/coverage/[system_short_name]-tumor/.", true,"auto");
$parser->addString("cov_pairs","ClinCNV file with tumor-normal pairs. Will be created from NGSD tumor-normal IDs by default.",true,"auto");
$parser->addInfile("system", "Processing system INI file (obligatory if NGSD is not available).", true);
$parser->addFlag("reanalyse_cohort","Reanalyse whole cohort of the same processing system.",true,"auto");
extract($parser->parse($argv));

//init
if(!db_is_enabled("NGSD") && ($n_id == "auto" || $cov_pairs == "auto" || !isset($system)))
{
	trigger_error("NGSD is not available. Please provide -n_id, -cov_pairs and -system manually.",E_USER_ERROR);
}

//creates file with all tumor-normal sample identifiers of the same processing system
function get_somatic_pairs($cov_folder_tumor,$file_name)
{
	$db = DB::getInstance("NGSD");
	$tumor_ids = glob($cov_folder_tumor ."/*.cov");
	for($i=0;$i<count($tumor_ids);++$i)
	{
		$tumor_ids[$i] = basename($tumor_ids[$i],".cov");
	}
	$tumor_normal_ids = array();
	
	foreach($tumor_ids as $tumor_id)
	{
		$t_sample_info = get_processed_sample_info($db,$tumor_id,false);
		if(is_null($t_sample_info)) continue; //skip samples that were not found in NGSD
		$normal_id = $t_sample_info["normal_name"];
		if(empty($normal_id)) continue;
		$tumor_normal_ids[$tumor_id] = $normal_id;
	}
	
	//write data to file
	$handle_tmp_file_pairs = fopen($file_name,'w');
	foreach($tumor_normal_ids as $tumor_id => $normal_id)
	{
		$line = "{$tumor_id},{$normal_id}\n";
		fwrite($handle_tmp_file_pairs,$line);
	}
	fclose($handle_tmp_file_pairs);
}

//sorts tsv file according to chromosome and start position
function sort_clincnv_file($in,$out)
{
	$data = file($in);
	$header = "";
	$contents ="";
	foreach($data as $line)
	{
		if(starts_with($line,"#"))
		{
			$header .= $line;
			continue;
		}
		$contents .= $line;
	}
	
	$tmp_file1 = temp_file(".tsv");
	$tmp_file2 = temp_file(".tsv");
	file_put_contents($tmp_file1,$contents);
	exec2("sort -k1.4,1 -k2,2 -V {$tmp_file1}  > {$tmp_file2}");
	
	$contents_sorted = file_get_contents($tmp_file2);
	
	$merged = $header . $contents_sorted;
	
	file_put_contents($out,$merged);
}

//Creates file with paths to coverage files using ref folder and current sample cov path
function create_file_with_paths($ref_cov_folder,$cov_path)
{
	$ref_paths = glob("{$ref_cov_folder}/*.cov");
	

	
	//Check whether cov file is already in ref cov folder -> remove from list
	for($i=0;$i<count($ref_paths);++$i)
	{
		if(strpos($ref_paths[$i],basename($cov_path,".cov")) !== false)
		{
			unset($ref_paths[$i]);
			$ref_paths = array_values($ref_paths); //reassign correct indices
			break;
		}
	}
	
	$ref_paths[] = $cov_path;

	for($i=0;$i<count($ref_paths);++$i)
	{
		$ref_paths[$i] .= "\n";
	}

	$out_file = temp_file(".txt");
	file_put_contents($out_file,$ref_paths);
	
	return $out_file;
}

if($n_id == "auto")
{
	$db = DB::getInstance("NGSD");
	$t_sample_info = get_processed_sample_info($db,$t_id,false);
	$n_id = $t_sample_info['normal_name'];
}

$data_folder = get_path("data_folder");
$sys = load_system($system,$t_id);

//get coverage folder (background)
if($cov_folder_n=="auto")
{
	$cov_folder_n = "{$data_folder}/coverage/".$sys['name_short'];
}
if(!is_dir($cov_folder_n))
{
	trigger_error("CNV calling skipped. Coverage files folder '$cov_folder_n' does not exist!", E_USER_WARNING);
	exit(0);
}
//get coverage tumor folder
if($cov_folder_t=="auto")
{
	$cov_folder_t = "{$data_folder}/coverage/".$sys['name_short']."-tumor";
}
if(!is_dir($cov_folder_t))
{
	trigger_error("CNV calling skipped. Coverage files folder '$cov_folder_t' does not exist!", E_USER_WARNING);
	exit(0);
}

$tmp_bed_annotated = $parser->tempFile(".bed");

$parser->exec(get_path("ngs-bits")."/"."BedAnnotateGC", " -in {$bed} -out {$tmp_bed_annotated}",true);
$parser->exec(get_path("ngs-bits")."/"."BedAnnotateGenes"," -in {$tmp_bed_annotated} -out {$tmp_bed_annotated}",true);


/************************
 * MERGE COVERAGE FILES *
 ************************/
$merge_files_exec_path = get_path("R_folder")."Rscript --vanilla ".get_path("ClinCnvSomatic")."/"."mergeFilesFromFolder.R";
$merged_cov_normal = $parser->tempFile(".txt");
$merged_cov_tumor = $parser->tempFile(".txt");

$cov_paths_n = create_file_with_paths($cov_folder_n,realpath($n_cov));
$cov_paths_t = create_file_with_paths($cov_folder_t,realpath($t_cov));

$parser->exec($merge_files_exec_path," -i {$cov_paths_n} -o {$merged_cov_normal}",true);
$parser->exec($merge_files_exec_path," -i {$cov_paths_t} -o {$merged_cov_tumor}",true);

/******************************************
 * CREATE PAIR FILE WITH TUMOR NORMAL IDS *
 ******************************************/
if($cov_pairs == "auto")
{
	$cov_pairs = $parser->tempFile(".txt");
	get_somatic_pairs($cov_folder_t,$cov_pairs);
}
//check whether coverage files of somatic pairs exist in cohort folder
$somatic_pairs = file($cov_pairs);
foreach($somatic_pairs as $line)
{
	if(starts_with($line,"#")) continue;
	if(empty(trim($line))) continue;
	list($pair_t,$pair_n) = explode(",",trim($line));
	
	if(strpos($pair_t,$t_id) !== false) continue; //current sample cov files can be stored in another folder than ref coverage files
	
	$cov_file_t = "{$cov_folder_t}/{$pair_t}.cov";
	$cov_file_n = "{$cov_folder_n}/{$pair_n}.cov" ;
	if(!file_exists($cov_file_t)) trigger_error("Could not find coverage file {$cov_file_t} for tumor ID {$pair_t}. Please check {$cov_pairs} file. Aborting.", E_USER_ERROR);
	if(!file_exists($cov_file_n)) trigger_error("Could not find coverage file {$cov_file_n} for normal ID {$pair_n}. Please check {$cov_pairs} file. Aborting.", E_USER_ERROR);
}

/*******************
 * EXECUTE CLINCNV *
 *******************/
//Determine folder for cohort output
if($cohort_folder == "auto")
{
	//create subdir in clincnv folder for cohort data
	if(!file_exists(get_path("ClinCnvSomatic") . "/cohorts"))
	{
		mkdir(get_path("ClinCnvSomatic") . "/cohorts");
		chmod($cohort_folder,0777);
	}
	$cohort_folder = get_path("ClinCnvSomatic") . "/" . "cohorts/".$sys['name_short'] . "/";
	if(!file_exists($cohort_folder))
	{
		mkdir($cohort_folder);
		chmod($cohort_folder,0777);
	}
}
if(!file_exists($cohort_folder))
{
	trigger_error("Directory for cohort output {$cohort_folder} does not exist.",E_USER_ERROR);
}
$call_cnvs_exec_path = get_path("R_folder")."Rscript --vanilla ".get_path("ClinCnvSomatic")."/"."firstStep.R";
$call_cnvs_params = " --normal {$merged_cov_normal} --tumor {$merged_cov_tumor} --out {$cohort_folder} --pair {$cov_pairs} --bed {$tmp_bed_annotated} --colNum 4 --folderWithScript ".get_path("ClinCnvSomatic")." --lengthS 1 --scoreS 40";
print($call_cnvs_params ."\n");
if($reanalyse_cohort) //Delete all sample files in cohort folder if reanalysis shall be performed
{
	if(is_dir("{$cohort_folder}/somatic/")) exec2("rm -r {$cohort_folder}/somatic/");
	if(is_dir("{$cohort_folder}/normal/")) exec2("rm -r {$cohort_folder}/normal/");
	$call_cnvs_params .=" --reanalyseCohort TRUE";
}
else
{
	$call_cnvs_params .=" --reanalyseCohort FALSE"; 
}

//Remove old data from somatic cohort folders
$cohort_folder_somatic_sample =  "{$cohort_folder}/somatic/{$t_id}-{$n_id}/";
if(is_dir($cohort_folder_somatic_sample)) exec2("rm -r {$cohort_folder_somatic_sample}");
$cohort_folder_normal_sample = "{$cohort_folder}/normal/{$n_id}/";
if(is_dir($cohort_folder_normal_sample)) exec2("rm -r {$cohort_folder_normal_sample}");


//Make sure all files in cohort folder have rw permisssions before exec
exec2("find {$cohort_folder} -type f -print0  | xargs -0 chmod 766");
$parser->exec($call_cnvs_exec_path,$call_cnvs_params,true);
exec2("find {$cohort_folder} -type f -print0  | xargs -0 chmod 766");


//copy segmentation files to folder containing output file
$cnvs_seg_file = "{$cohort_folder_somatic_sample}/{$t_id}-{$n_id}_cnvs.seg";
if(file_exists($cnvs_seg_file))
{
	$parser->copyFile($cnvs_seg_file,dirname($out)."/{$t_id}-{$n_id}_cnvs.seg");
}
$cnvs_cov_file = "{$cohort_folder_somatic_sample}/{$t_id}-{$n_id}_cov.seg";
if(file_exists($cnvs_cov_file))
{
	$parser->copyFile($cnvs_cov_file,dirname($out)."/{$t_id}-{$n_id}_cov.seg");
}

//check whether desired tumor normal pair was created successfully
$unparsed_cnv_file = "{$cohort_folder_somatic_sample}/CNAs.txt";
if(!file_exists($unparsed_cnv_file))
{
	trigger_error("ClinCNV output file {$unparsed_cnv_file} was not created. Aborting.",E_USER_ERROR);
}

/**************************
 * PARSE RAW CLINCNV FILE *
 **************************/
//sort by chromosome and cnv start position
sort_clincnv_file($unparsed_cnv_file,$unparsed_cnv_file);
 
//insert column with sample identifier for convenience (next to "end"-column")
$cnvs = Matrix::fromTSV($unparsed_cnv_file);

$sample_col = array_fill(0,$cnvs->rows(),"{$t_id}-{$n_id}");
$cnvs->insertCol($cnvs->getColumnIndex("end")+1,$sample_col,"sample");

//insert column with CNV sizes
$i_start = $cnvs->getColumnIndex("start");
$i_end = $cnvs->getColumnIndex("end");
$col_sizes = array();
for($r=0;$r<$cnvs->rows();++$r)
{
	$start = $cnvs->get($r,$i_start);
	$end = $cnvs->get($r,$i_end);
	$size = $end - $start + 1;
	$col_sizes[] = $size;
}

$cnvs->insertCol($cnvs->getColumnIndex("sample")+1,$col_sizes,"size");

/**********************
 * DETERMINE CNV TYPE *
 **********************/
$chr = "";
$start = 0;
$end = 0;
$arm = "";
$chrom_arms = array();
$cyto_bands = Matrix::fromTSV(repository_basedir()."/data/dbs/UCSC/cytoBand.txt");

for($i=0;$i<$cyto_bands->rows();++$i)
{
	$row = $cyto_bands->getRow($i);
	list($tmp_chr,$tmp_start,$tmp_end,$tmp_arm) = $row;
	++$tmp_start; // 0 to 1-base
	$tmp_arm = $tmp_arm[0]; //p-arm or q-arm	
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
$idx_genes = $cnvs->getColumnIndex("genes");
for($i=0;$i<$cnvs->rows();++$i)
{
	$row = $cnvs->getRow($i);
	
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
$cnvs->addCol($new_col, "cnv_type","Type of CNV: focal (< 25 % of chrom. arm, < 3 genes), cluster (focal, > 3 genes), (partial) p/q-arm, (partial) chromosome.");
$cnvs->toTSV($out);
?>