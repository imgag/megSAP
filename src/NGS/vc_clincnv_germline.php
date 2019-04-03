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
$parser->addInt("cov_max", "Maximum number of coverage files used for CNV analysis, used to keep run-time and RAM requirement manageable (~TODO for WES and ~TODO for WGS).", true, 200);
$parser->addInt("max_cnvs", "Number of expected CNVs (~200 for WES and ~2000 for WGS).", true, 2000);

extract($parser->parse($argv));

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
	//$t_start = microtime(true);
	$output = array();

	$i = -1;
	$h = fopen($filename, "r");
	while(!feof($h))
	{
		$line = fgets($h);
		
		++$i;
		if ($rows_to_use[$i]==false) continue;
		
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		
		list($chr, $start, $end, $avg_cov) = explode("\t", $line);
		$output[$chr][] = $avg_cov;
	}
	//print "COV DONE ".time_readable(microtime(true) - $t_start)." (".basename($filename).")\n";
}

//init
$ps_name = basename($cov,".cov");

//determine coverage files
$cov_files = glob($cov_folder."/*.cov");
$cov_files[] = $cov;
$cov_files = array_unique(array_map("realpath", $cov_files));
if (count($cov_files)<$cov_min)
{
	trigger_error("CNV calling skipped. Only ".count($cov_files)." coverage files found in folder '$cov_folder'. At least {$cov_min} files are needed!", E_USER_WARNING);
	exit(0);
}

//select the most similar samples
if (count($cov_files)>$cov_max)
{
	$parser->log(count($cov_files)." coverge files are given. Selecting the most similar {$cov_max} samples as reference samples!");
	
	//sort coverage files
	sort($cov_files);
	
	//create target region without polymorphic regions
	//$t_start = microtime(true);
	$poly_merged = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedAdd", "-in ".repository_basedir()."/data/misc/cn_polymorphisms_zarrei.bed ".repository_basedir()."/data/misc/cn_polymorphisms_demidov.bed ".repository_basedir()."/data/misc/centromer_telomer_hg19.bed -out {$poly_merged}", true);
	$roi_poly = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedIntersect", "-in {$bed} -in2 {$poly_merged} -out {$roi_poly} -mode in", true);
	$roi_nonpoly = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedSubtract", "-in {$bed} -in2 {$roi_poly} -out {$roi_nonpoly}", true);
	//print "PREP DONE ".time_readable(microtime(true) - $t_start)."\n";

	//determine which rows of coverage profiles to use
	//$t_start = microtime(true);
	$rows_to_use = determine_rows_to_use($cov, $roi_nonpoly);
	//print "INDICES DONE ".time_readable(microtime(true) - $t_start)."\n";
	
	//determine correlation of each sample/cov file
	$file2corr = array();
	$cov1 = null;
	$cov2 = null;
	load_coverage_profile($cov, $rows_to_use, $cov1);
	foreach($cov_files as $cov_file)
	{
		load_coverage_profile($cov_file, $rows_to_use, $cov2);
		
		//$t_start = microtime(true);
		$corr = array();
		foreach($cov1 as $chr => $profile1)
		{
			$corr[] = number_format(correlation($profile1, $cov2[$chr]), 2);
		}
		//print "CALC DONE ".time_readable(microtime(true) - $t_start)."\n";
		
		$file2corr[$cov_file] = number_format(median($corr),2);
	}
		
	//sort by correlation
	arsort($file2corr);
	
	//selects best n
	$file2corr = array_slice($file2corr, 0, $cov_max);
	$add_info = array();
	foreach ($file2corr as $f => $c)
	{
		$add_info[] = "$f ($c)";
	}
	$parser->log("Selected the following files as reference samples (correlation):", $add_info);
	
	$cov_files = array_keys($file2corr);
}

//merge coverage files to one file
$prefix = "";
$max_open_files_soft = trim(implode("", exec2("ulimit -Sn")[0]));
if (count($cov_files)>=$max_open_files_soft)
{
	$max_open_files_hard = trim(implode("", exec2("ulimit -Hn")[0]));
	if (count($cov_files)>=$max_open_files_hard)
	{
		trigger_error("Number of coverage files (".count($cov_files).") is bigger than the maximum number of open files (hard limit: {$max_open_files_hard})!", E_USER_ERROR);
	}
	else
	{
		trigger_error("Number of coverage files ({".count($cov_files)."}) is bigger than the maximum number of open files (soft limit: {$max_open_files_soft}). Trying to increase the soft limit...", E_USER_WARNING);
		$prefix = "ulimit -Sn 10000; ";
	}
}
$tmp = $parser->tempFile(".txt");
file_put_contents($tmp, implode("\n", $cov_files));
$cov_merged = $parser->tempFile(".cov");
$parser->exec($prefix.get_path("ngs-bits")."TsvMerge", "-in $tmp -cols chr,start,end -simple -out {$cov_merged}", true);

//execute ClinCNV
$out_folder = $parser->tempFolder();
$args = [
"--normal {$cov_merged}",
"--normalSample {$ps_name}",
"--bed {$bed}",
"--out {$out_folder}",
"--maxNumGermCNVs {$max_cnvs}",
"--numberOfThreads {$threads}",
"--lengthG 0", //actually means length 1 ;)
"--scoreG 20",
"--superRecall 3"
];
$parser->exec(get_path("clincnv")."/clinCNV.R", implode(" ", $args), true);


//sort and extract sample data from output folder
$parser->exec(get_path("ngs-bits")."/BedSort","-in {$out_folder}/normal/{$ps_name}/{$ps_name}_cnvs.tsv -out $out",true);
$parser->copyFile("{$out_folder}/normal/{$ps_name}/{$ps_name}_cov.seg", substr($out, 0, -4).".seg");
$parser->copyFile("{$out_folder}/normal/{$ps_name}/{$ps_name}_cnvs.seg", substr($out, 0, -4)."_cnvs.seg");

//Add header line holding type of analysis
file_put_contents($out,"##ANALYSISTYPE=CLINCNV_GERMLINE_SINGLE\n".file_get_contents($out));

//annotate
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 ".repository_basedir()."/data/misc/cn_polymorphisms_zarrei.bed -out {$out}", true);
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 ".repository_basedir()."/data/gene_lists/dosage_sensitive_disease_genes.bed -out {$out}", true);
$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; 
if (file_exists($omim_file)) //optional because of license
{
	$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$omim_file} -out {$out}", true);
}

?>

