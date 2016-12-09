<?php

/**
	@page x_diagnostics
*/

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";

require_once($basedir."Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("x_diagnostics", "\$Rev: 849 $", "XLMR diagnostics pipeline.");
$parser->addInfile("bam",  "X chromosome BAM file.", false);
$parser->addInfile("system",  "Processing system INI file.", false);
$parser->addString("out_folder", "Output folder.", false);
extract($parser->parse($argv));

//system
$sys = load_system($system);

//check output directory
if(!is_dir($out_folder) || !is_writeable($out_folder))
{
	trigger_error("Cannot access output folder '$out_folder'!", E_USER_ERROR);
}
$basename = basename($bam, ".bam");

// determine gender
list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $bam", true);
$gender = trim(substr($stdout[count($stdout)-1], 7));
$rpt = array();
$rpt[] = "Geschlecht: ".($gender=="female" ? "w" : "m");

// detect exon deletions
if ($gender!="female")
{
	//create cov file
	$cov_folder = get_path("data_folder")."/coverage/".$sys['name_short']."/";
	$cov_file = $cov_folder.$basename.".bam.cov";
	$parser->exec(get_path("ngs-bits")."BedCoverage", "-in ".$sys["target_file"]." -bam $bam -out $cov_file -min_mapq 0", true);
	chmod($cov_file, 0775);

	//run CnvHunter
	$temp1 = $parser->tempFile(".tsv");
	list($stdout) = $parser->exec(get_path("ngs-bits")."CnvHunter", "-in ".$cov_folder."*.cov -out $temp1 -sam_min_corr 0.85", true);
	
	//check if the current sample is invalid
	$bad_sample = "";
	foreach($stdout as $line)
	{
		if (contains($line, "bad sample:") && contains($line, $basename))
		{
			$bad_sample = trim(strtr($line, array("bad sample:"=>"")));
			break;
		}
	}
		
	//create CNV output
	$rpt[] = "";
	$rpt[] = "CNV-Detection:";
	if ($bad_sample!="")
	{
		$rpt[] = "  Nicht möglich weil: $bad_sample";
	}
	else
	{
		$cnvs = array();
		$file = load_tsv($temp1);
		foreach($file as $line)
		{
			if ($line[3]==$basename)
			{
				$cnvs[] = $line;
			}
		}
		$rpt[] = "  Mögliche CNVs: ".count($cnvs);
		$rpt[] = "  #Koordinaten\tExons\tCopyNumber\tDetails";
		foreach($cnvs as $cnv)
		{
			$rpt[] = "  ".$cnv[0].":".$cnv[1]."-".$cnv[2]."\t".$cnv[4]."\t".$cnv[5]."\t".$cnv[6];
		}
	}
}

//write basic report
file_put_contents($out_folder."/{$basename}_".$sys['name_short']."_report.txt", implode("\n", $rpt));

//XLMR low-coverage regions
if ($sys['name_short']=="ssX")
{
	$roi = get_path("data_folder")."/gene_lists/XLMR_known_genes_cds_merged_v4.bed";
	$cutoff = ($gender=="female") ? 20 : 10;
	$lowcov_file = $out_folder."/{$basename}_XLMR_lowcov.bed";
	$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $roi -bam $bam -out $lowcov_file -cutoff $cutoff", false);
	$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $lowcov_file -extend 25 -out $lowcov_file", true);
}

?>
