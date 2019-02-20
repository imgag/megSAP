<?php 
/** 
	@page vc_clincnv_germline
	
	@todo update copy_number_map_strict
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_clincnv_germline", "Germline copy-number calling with ClinCNV.");
$parser->addInfile("cov", "Coverage file for sample (tab-separated file with columns chr, start, end, coverage).", false);
$parser->addInfile("cov_folder", "Coverage files folder.", false);
$parser->addInfile("bed", "BED file with annotations e.g. GC-content and gene names.", false);
$parser->addOutFile("out", "Output file in TSV format.", false);
$parser->addInt("max_cnvs", "Number of expected CNVs (~100 for WES and ~TODO for WGS).", false);
//optional
$parser->addInt("cov_min", "Minimum number of coverage files required for CNV analysis.", true, 20);

extract($parser->parse($argv));

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

//merge coverage files to one file
$tmp = $parser->tempFile(".txt");
file_put_contents($tmp, implode("\n", $cov_files));
$cov_merged = $parser->tempFile(".cov");
$parser->exec(get_path("ngs-bits")."TsvMerge", "-in $tmp -cols chr,start,end -simple -out {$cov_merged}", true);

//execute ClinCNV
$out_folder = $parser->tempFolder();
$args = [
"--normal {$cov_merged}",
"--normalSample {$ps_name}",
"--bed {$bed}",
"--out {$out_folder}",
"--maxNumGermCNVs {$max_cnvs}"
];
$parser->exec(get_path("clincnv")."/clinCNV.R", implode(" ", $args), true);

//sort and extract sample data from output folder
$parser->exec(get_path("ngs-bits")."/BedSort","-in {$out_folder}/normal/{$ps_name}/{$ps_name}_cnvs.tsv -out $out",true);
$parser->copyFile("{$out_folder}/normal/{$ps_name}/{$ps_name}_cov.seg", substr($out, 0, -4).".seg");

//annotate
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 ".repository_basedir()."/data/misc/copy_number_map_strict.bed -out {$out}", true);
$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 ".repository_basedir()."/data/gene_lists/dosage_sensitive_disease_genes.bed -out {$out}", true);
$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; 
if (file_exists($omim_file)) //optional because of license
{
	$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$out} -in2 {$omim_file} -out {$out}", true);
}

?>

