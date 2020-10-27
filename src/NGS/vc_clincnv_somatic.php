<?php 
/** 
	@page vc_clincnv_somatic

*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_clincnv_somatic", "Wrapper for CnvHunter tool.");
$parser->addString("t_id", "Processed sample id of tumor, e.g. 'GS123456_01'.", false);
$parser->addInfile("t_cov","Coverage file for tumor sample",false);
$parser->addInfile("n_cov","Coverage file for normal sample",false);
$parser->addOutFile("out", "Output file.",false);
//optional
$parser->addInfile("bed","Bed file for target region.",false);
$parser->addInfile("bed_off","Off-target bed file.",true); //s_dna
$parser->addInfile("t_cov_off","Off-target coverage file for tumor sample",true);
$parser->addInfile("n_cov_off","Off-target coverage file for normal sample",true);
$parser->addString("cov_folder_n_off", "Folder with off-target normal coverage files (if different from [data_folder]/coverage/[system_short_name]_off_target/.", true, "auto");
$parser->addString("cov_folder_t_off", "Folder with all tumor coverage files (if different from [data_folder]/coverage/[system_short_name]-tumor_off_target/.", true,"auto");
$parser->addString("cohort_folder", "Folder that contains ClinCNV cohort data (same processing system).", true,"auto");
$parser->addString("n_id","Processed sample id of normal sample, e.g. 'GS123456_01'. Using normal sample ID from NGSD by default.",true,"auto");
$parser->addString("cov_folder_n", "Folder with all coverage files (if different from [data_folder]/coverage/[system_short_name]/.", true, "auto");
$parser->addString("cov_folder_t", "Folder with all tumor coverage files (if different from [data_folder]/coverage/[system_short_name]-tumor/.", true,"auto");
$parser->addString("cov_pairs","ClinCNV file with tumor-normal pairs. Will be created from NGSD tumor-normal IDs by default.",true,"auto");
$parser->addString("baf_folder","Folder containing files with B-Allele frequencies.",true);
$parser->addInfile("system", "Processing system INI file (obligatory if NGSD is not available).", true);
$parser->addFlag("reanalyse_cohort","Reanalyse whole cohort of the same processing system.",true,"auto");
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("guide_baseline","baseline region, format e.g. chr1:12-12532",true);
extract($parser->parse($argv));

//init
if(!db_is_enabled("NGSD") && ($n_id == "auto" || $cov_pairs == "auto" || !isset($system)))
{
	trigger_error("NGSD is not available. Please provide -n_id, -cov_pairs and -system manually.", E_USER_ERROR);
}

$use_off_target = true;
if(!isset($bed_off) || !isset($t_cov_off) || !isset($n_cov_off))
{
	$use_off_target = false;
}

//extract tool path from ClinCNV command
$tool_folder = "";
$parts = explode(" ", get_path("clincnv"));
foreach($parts as $part)
{
	if ($part!="" && file_exists($part))
	{
		$tool_folder = $part;
	}
}
if ($tool_folder=="")
{
	trigger_error("Could not determine tool folder of ClinCNV!", E_USER_ERROR);
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
	$handle_tmp_file_pairs = fopen2($file_name,'w');
	foreach($tumor_normal_ids as $tumor_id => $normal_id)
	{
		$line = "{$tumor_id},{$normal_id}\n";
		fwrite($handle_tmp_file_pairs,$line);
	}
	fclose($handle_tmp_file_pairs);
}

//Creates file with paths to coverage files using ref folder and current sample cov path, if $sample_ids is set: skip all ids not contained
function create_file_with_paths($ref_cov_folder,$cov_path, $sample_ids = array())
{
	global $parser;
	
	$ref_paths = glob("{$ref_cov_folder}/*.cov");
	
	$paths_to_be_included = array();
	

	//Remove samples that are not in $sample_ids
	if(count($sample_ids) > 0)
	{	
		for($i=0;$i<count($ref_paths);++$i)
		{	
			$id = basename($ref_paths[$i], ".cov");

			if(in_array($id, $sample_ids))
			{
				$paths_to_be_included[] = $ref_paths[$i];
			}
		}
	}
	else
	{
		for($i=0;$i<count($ref_paths);++$i)
		{	
			$paths_to_be_included[] = $ref_paths[$i];
		}
	}

	//Check whether cov file is already in cov folder -> remove from list (could be specified in another dir!)
	for($i=0;$i<count($paths_to_be_included);++$i)
	{
		if(strpos($paths_to_be_included[$i],basename($cov_path,".cov")) !== false)
		{
			unset($paths_to_be_included[$i]);
			$paths_to_be_included = array_values($paths_to_be_included); //reassign correct indices
			break;
		}
	}
	$paths_to_be_included[] = $cov_path;
	

	$out_file = $parser->tempFile(".txt");
	file_put_contents($out_file,implode("\n",$paths_to_be_included) );
	
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
	trigger_error("CNV calling skipped. Coverage files folder cov_folder_n '$cov_folder_n' does not exist!", E_USER_WARNING);
	exit(0);
}
if($cov_folder_n_off == "auto")
{
	$cov_folder_n_off = "{$data_folder}/coverage/".$sys['name_short']."_off_target";
}
if($use_off_target && !is_dir($cov_folder_n_off))
{
	trigger_error("CNV calling skipped. Off-target Coverage files folder cov_folder_n_off '$cov_folder_n_off' does not exist!", E_USER_WARNING);
	exit(0);
}

//get coverage tumor folder
if($cov_folder_t=="auto")
{
	$cov_folder_t = "{$data_folder}/coverage/".$sys['name_short']."-tumor";
}
if(!is_dir($cov_folder_t))
{
	trigger_error("CNV calling skipped. Coverage files folder cov_folder_t '$cov_folder_t' does not exist!", E_USER_WARNING);
	exit(0);
}
if($cov_folder_t_off=="auto")
{
	$cov_folder_t_off = "{$data_folder}/coverage/".$sys['name_short']."-tumor_off_target";
}
if($use_off_target && !is_dir($cov_folder_t_off))
{
	trigger_error("CNV calling skipped. Coverage files folder cov_folder_t_off '$cov_folder_t_off' does not exist!", E_USER_WARNING);
	exit(0);
}

$tmp_bed_annotated = $parser->tempFile(".bed");

$parser->exec(get_path("ngs-bits")."/"."BedAnnotateGC", " -in {$bed} -out {$tmp_bed_annotated} -ref ".genome_fasta($sys['build']),true);
$parser->exec(get_path("ngs-bits")."/"."BedAnnotateGenes"," -in {$tmp_bed_annotated} -out {$tmp_bed_annotated}",true);


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

//Make list of tumor ids and normal ids stored in tsv list 
$sample_tids = array();
$sample_nids = array();
foreach($somatic_pairs as $line)
{
	if(starts_with($line,"#")) continue;
	if(empty(trim($line))) continue;
	list($tid,$nid) = explode(",",trim($line));
	$sample_tids[] = $tid;
	$sample_nids[] = $nid;
}

/************************
 * MERGE COVERAGE FILES *
 ************************/
$merged_cov_normal = $parser->tempFile(".txt");
$merged_cov_tumor = $parser->tempFile(".txt");
$cov_paths_n = create_file_with_paths($cov_folder_n,realpath($n_cov), $sample_nids);
$cov_paths_t = create_file_with_paths($cov_folder_t,realpath($t_cov), $sample_tids);

$parser->exec(get_path("ngs-bits")."/TsvMerge" , " -in $cov_paths_n -out {$merged_cov_normal} -cols chr,start,end -simple",true);
$parser->exec(get_path("ngs-bits")."/TsvMerge" , " -in $cov_paths_t -out {$merged_cov_tumor} -cols chr,start,end -simple",true);

//off-targets
if($use_off_target)
{
	$merged_cov_tumor_off = $parser->tempFile(".txt");
	$merged_cov_normal_off = $parser->tempFile(".txt");
	$cov_paths_n_off = create_file_with_paths($cov_folder_n_off,realpath($n_cov_off), $sample_nids);
	$cov_paths_t_off = create_file_with_paths($cov_folder_t_off,realpath($t_cov_off), $sample_tids);
	$parser->exec(get_path("ngs-bits")."/TsvMerge" , " -in $cov_paths_n_off -out {$merged_cov_normal_off} -cols chr,start,end -simple",true, $sample_nids);
	$parser->exec(get_path("ngs-bits")."/TsvMerge" , " -in $cov_paths_t_off -out {$merged_cov_tumor_off} -cols chr,start,end -simple",true, $sample_tids);
}



/*******************
 * EXECUTE CLINCNV *
 *******************/
//Determine folder for cohort output, standard is in ClinCNV app directory.
if($cohort_folder == "auto")
{	
	//create main cohorts folder
	if(!file_exists("{$tool_folder}/cohorts"))
	{
		mkdir("{$tool_folder}/cohorts");
		chmod("{$tool_folder}/cohorts", 0777);
	}
	
	//create system-specific cohorts folder
	$cohort_folder = "{$tool_folder}/cohorts/".$sys['name_short']."/";
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
$call_cnvs_exec_path = get_path("clincnv")."/clinCNV.R";

$args = [
"--normal", $merged_cov_normal,
"--tumor", $merged_cov_tumor,
"--out", $cohort_folder,
"--pair", $cov_pairs,
"--bed", $tmp_bed_annotated,
"--colNum", "4",
"--lengthS", "5",
"--scoreS", "100",
"--numberOfThreads {$threads}"
];

if($use_off_target)
{
	$args[] = "--tumorOfftarget $merged_cov_tumor_off";
	$args[] = "--normalOfftarget $merged_cov_normal_off";
	$args[] = "--bedOfftarget $bed_off";
}

if(isset($guide_baseline))
{
	$args[] = "--guideBaseline $guide_baseline";
}

if($reanalyse_cohort) //Delete all sample files in cohort folder if reanalysis shall be performed
{
	if(is_dir("{$cohort_folder}/somatic/")) exec2("rm -r {$cohort_folder}/somatic/");
	if(is_dir("{$cohort_folder}/normal/")) exec2("rm -r {$cohort_folder}/normal/");
	$args[] = "--reanalyseCohort TRUE";
}
else //specify single output sample otherwise
{
	$args[] = "--reanalyseCohort FALSE"; 
	$args[] = "--normalSample {$n_id}";
	$args[] = "--tumorSample {$t_id}";
}

if(is_dir($baf_folder)) $args[] = "--bafFolder {$baf_folder}";

//Remove old data from somatic cohort folders
$cohort_folder_somatic_sample =  "{$cohort_folder}/somatic/{$t_id}-{$n_id}/";
if(is_dir($cohort_folder_somatic_sample)) exec2("rm -r {$cohort_folder_somatic_sample}");
$cohort_folder_normal_sample = "{$cohort_folder}/normal/{$n_id}/";
if(is_dir($cohort_folder_normal_sample)) exec2("rm -r {$cohort_folder_normal_sample}");

//execite ClinCNV and make sure all files in cohort folder have 777 permisssions
function chmod_recursive($folder)
{
	list(, , $exit_code) = exec2("chmod -R 777 {$folder}", false);
	if ($exit_code!=0)
	{
		trigger_error("Could not change permissions of all files in {$folder}!", E_USER_WARNING);
	}
}
chmod_recursive($cohort_folder);

$parser->exec($call_cnvs_exec_path,implode(" ",$args),true);
chmod_recursive($cohort_folder);

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

$cnv_plots = glob("{$cohort_folder_somatic_sample}/*.png");
if(count($cnv_plots) > 0)
{
	foreach($cnv_plots as $plot)
	{
		$parser->copyFile($plot,dirname($out). "/" .basename($plot));
	}
}

//check whether desired tumor normal pair was created successfully
$unparsed_cnv_file = "{$cohort_folder_somatic_sample}/CNAs_{$t_id}-{$n_id}.txt";
if(!file_exists($unparsed_cnv_file))
{
	trigger_error("ClinCNV output file {$unparsed_cnv_file} was not created. Aborting.",E_USER_ERROR);
}

/**************************
 * PARSE RAW CLINCNV FILE *
 **************************/
//sort by chromosome and cnv start position
$parser->exec(get_path("ngs-bits")."/BedSort","-in $unparsed_cnv_file -out $unparsed_cnv_file",true);
 
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

//remove unneeded columns
$cnvs->removeCol($cnvs->getColumnIndex("median_loglikelihood"));
$cnvs->removeCol($cnvs->getColumnIndex("major_CN_allele2"));
$cnvs->removeCol($cnvs->getColumnIndex("minor_CN_allele2"));
$cnvs->removeCol($cnvs->getColumnIndex("tumor_clonality2"));


//Add baseline value as comment to output file
if(isset($guide_baseline))
{
	$cnvs->addComment("#guideBaseline: $guide_baseline");
}


/**********************
 * DETERMINE CNV TYPE *
 **********************/
$chr = "";
$start = 0;
$end = 0;
$arm = "";
$chrom_arms = array();
$cyto_bands = Matrix::fromTSV(repository_basedir()."/data/misc/cytoBand.txt");

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


//store output tsv file
$cnvs->toTSV($out);

//annotate additional gene info
$parser->exec(get_path("ngs-bits")."CnvGeneAnnotation", "-in {$out} -out {$out}", true);
//annotate overlap with pathogenic CNVs
if(db_is_enabled("NGSD")) $parser->exec(get_path("ngs-bits")."NGSDAnnotateCNV", "-in {$out} -out {$out}", true);

?>