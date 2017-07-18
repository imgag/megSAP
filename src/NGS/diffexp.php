<?php

/**
 @page diffexp
 */
require_once (dirname($_SERVER['SCRIPT_FILENAME']) . "/../Common/all.php");

$parser = new ToolBase("diffexp", "Run a differential expression analysis.");

$parser->addStringArray("samples_case", "Samples in the case group.", false);
$parser->addStringArray("samples_ctrl", "Samples in the control group.", false);

$parser->addString("case_name", "Name of the case group.", true, "case");
$parser->addString("ctrl_name", "Name of the control group.", true, "ctrl");

$parser->addFloat("filter_pval", "Cut-off for the p-value filter.", true, 0.05);
$parser->addFloat("filter_fdr", "Cut-off for the FDR filter.", true, 0.05);
$parser->addFloat("filter_fc", "Cut-off for the log fold change filter.", true, 1);

$parser->addFlag("prepare_targets", "Only prepare targets file for analysis script but do not run the analysis.");

extract($parser->parse($argv));

// targets file
$targets = new Matrix();
$colnames = array(
	"Sample",
	"Name",
	"Group",
	"File"
);
// $targets->setHeaders($colnames);
$targets->addRow($colnames);

// check existance of each sample
// TODO check that all samples have same processing system/build (?)
foreach (array_merge($samples_case, $samples_ctrl) as $sample)
{
	$info = get_processed_sample_info($sample);
	$sample_folder = $info["ps_folder"];
	
	// check that raw read counts files are there
	$readcounts = $sample_folder . $sample . "_counts_raw.tsv";
	if (! is_file($readcounts))
	{
		trigger_error("Could not find raw read counts file for sample '$sample'!", E_USER_ERROR);
	}
	
	$group = in_array($sample, $samples_case) ? $case_name : $ctrl_name;
	$targets->addRow(array(
		$sample,
		$info["name_external"],
		$group,
		$readcounts
	));
}

// contrasts file
$contrasts = new Matrix();
$contrasts->addRow(array(
	"Case",
	"Control",
	"Name"
));
$contrasts->addRow(array(
	$case_name,
	$ctrl_name,
	""
));

$temp = $parser->tempFolder("diffexp");
$targets_f = $temp . "/targets.tsv";
$targets->toTSV($targets_f);

if ($prepare_targets)
{
	print(file_get_contents($targets_f));
}
else
{
	$build = "TODO";
	$contrasts_f = $temp . "/contrasts.tsv";
	$contrasts->toTSV($contrasts_f);
	$args = array(
		"Project",
		$targets_f,
		$contrasts_f,
		"hsapiens_gene_ensembl",
		"ensembl_gene_id"
	);
	$parser->exec("Rscript ~/ahmattj1/Rcode/expression/expr_cli.R", implode(" ", $args), true);
}
