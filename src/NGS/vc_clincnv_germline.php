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
//optional
$parser->addInt("cov_min", "Minimum number of coverage files required for CNV analysis.", true, 30);
extract($parser->parse($argv));

//init
$ps_name = basename($cov,".cov");

//extract tool path from ClinCNV command
$tool_folder = "";
$parts = explode(" ", get_path("clincnv_germline"));
foreach($parts as $part)
{
	if ($part!="" && file_exists($part))
	{
		$tool_folder = $part;
	}
}
if ($tool_folder=="") trigger_error("Could not determine tool folder of ClinCNV!", E_USER_ERROR);

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
"--folderWithScript {$tool_folder}"
];
$parser->exec(get_path("clincnv_germline")."/firstStep.R", implode(" ", $args), true);

//extract sample data from output folder
$parser->copyFile("{$out_folder}/normal/{$ps_name}/{$ps_name}_cnvs.tsv", $out);
$parser->copyFile("{$out_folder}/normal/{$ps_name}/{$ps_name}_cov.seg", substr($out, 0, -4).".seg");

//TODO:
//- annotate CNP regions, dosage-sensitive genes, OMIM
//- use same version in germline and somatic

?>

