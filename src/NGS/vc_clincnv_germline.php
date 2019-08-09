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
$parser->addInt("max_tries", "Maximum number of tries for calling ClinCNV (R parallelization sometimes breaks with no reason", true, 10);

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
	$h = fopen2($filename, "r");
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
$repository_basedir = repository_basedir();
$data_folder = get_path("data_folder");
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
	$parser->exec(get_path("ngs-bits")."BedAdd", "-in {$repository_basedir}/data/misc/cnps_genomes_imgag.bed {$repository_basedir}/data/misc/centromer_telomer_hg19.bed -out {$poly_merged}", true);
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
$tmp = $parser->tempFile(".txt");
file_put_contents($tmp, implode("\n", $cov_files));
$cov_merged = $parser->tempFile(".cov");
$parser->exec(get_path("ngs-bits")."TsvMerge", "-in $tmp -cols chr,start,end -simple -out {$cov_merged}", true);

//execute ClinCNV (with workaround for hanging jobs)
$out_folder = $parser->tempFolder();
$args = [
"--normal {$cov_merged}",
"--normalSample {$ps_name}",
"--bed {$bed}",
"--out {$out_folder}",
"--maxNumGermCNVs {$max_cnvs}",
"--numberOfThreads {$threads}",
"--lengthG 1", //actually means length 2 - superRecall will call down to one region in the second step 
"--scoreG 20",
"--superRecall 3"
];

$command = get_path("clincnv")."/clinCNV.R";
$parameters = implode(" ", $args);
$pid = getmypid();
$stdout_file = $parser->tempFile(".stdout", "megSAP_clincnv_pid{$pid}_");
$stderr_file = $parser->tempFile(".stderr", "megSAP_clincnv_pid{$pid}_");
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

//sort and extract sample data from output folder
$parser->exec(get_path("ngs-bits")."/BedSort","-in {$out_folder}/normal/{$ps_name}/{$ps_name}_cnvs.tsv -out $out",true);
$parser->copyFile("{$out_folder}/normal/{$ps_name}/{$ps_name}_cov.seg", substr($out, 0, -4).".seg");
$parser->copyFile("{$out_folder}/normal/{$ps_name}/{$ps_name}_cnvs.seg", substr($out, 0, -4)."_cnvs.seg");

//add high-quality CNV count to header
$cnv_calls = file($out);
$hq_cnvs = 0;
foreach($cnv_calls as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	list($c, $s, $e, $tmp_cn, $ll) = explode("\t", $line);
	if ($ll>=20)
	{
		++$hq_cnvs;
	}	
}
array_splice($cnv_calls, 3, 0, array("##high-quality cnvs: {$hq_cnvs}\n"));
file_put_contents($out, $cnv_calls);

//annotate
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$repository_basedir}/data/misc/cnps_genomes_imgag.bed -overlap -out {$out}", true);
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$repository_basedir}/data/misc/cnps_700genomes_pcawg.bed -overlap -out {$out}", true);
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -out {$out}", true);
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes.bed -no_duplicates -out {$out}", true);
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs.bed -no_duplicates -out {$out}", true);
$hgmd_file = "{$data_folder}/dbs/HGMD/hgmd_cnvs.bed"; //optional because of license
if (file_exists($hgmd_file))
{
	$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$hgmd_file} -no_duplicates -out {$out}", true);
}
$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
if (file_exists($omim_file))
{
	$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$omim_file} -no_duplicates -out {$out}", true);
}

?>

