<?php

/**
  @page vc_manta
  @TODO improve config file handling
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_manta", "Call somatic strucural variants with manta. Creates an VCF file.");
$parser->addInfileArray("bam", "Normal BAM file(s). Only one BAM file allowed for somatic mode.", false);
$parser->addString("out", "Output file (gzipped and tabix indexed).", false);
//optional
$parser->addInfile("t_bam", "Tumor BAM file for somatic mode.", true);
$parser->addFlag("exome", "If set, settings for exome are used.", true);
$parser->addString("smallIndels", "Output file for candidate small indels.", true, "");
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addString("temp", "Temporary folder.", true, "auto");
$parser->addString("regions", "Comma-separated string of regions.", true, "");

$parser->addEnum("config_preset", "Use preset configuration.", true, array("default", "high_sensitivity"), "default");

$parser->addInt("threads", "Number of threads used.", true, 4);

extract($parser->parse($argv));

//determine somatic mode
$somatic_mode = isset($t_bam);
if ($somatic_mode && count($bam) > 1) trigger_error("Only one bam file allowed in somatic mode!", E_USER_ERROR);

//resolve configuration preset
$config_default = get_path("manta")."/configManta.py.ini";
$config = $config_default;
if ($config_preset === "high_sensitivity")
{
	$config = $parser->tempFile();
	//change minEdgeObservations and minCandidateSpanningCount parameters to 2
	$parser->exec("sed", "'s/^minEdgeObservations.\\+$/minEdgeObservations = 2/; s/^minCandidateSpanningCount.\\+/minCandidateSpanningCount = 2/' {$config_default} > {$config}", false);
	$parser->log("Manta config file:", file($config));
}

$temp_folder = $temp === "auto" ? $parser->tempFolder() : $temp;
$manta_folder = "{$temp_folder}/mantaAnalysis";
$genome = get_path("local_data")."/{$build}.fa";

$pars = "";
$pars .= "--referenceFasta {$genome} ";
$pars .= "--runDir {$manta_folder} --config {$config} ";
$pars .= "--normalBam ".implode(" --normalBam ", $bam)." ";
if ($somatic_mode) $pars .= "--tumorBam $t_bam ";
if ($exome) $pars .= "--exome ";
if (!empty($regions)) $pars .= " --region ".implode(" --region ", explode(",", $regions))." ";

//run manta
$parser->exec("python ".get_path('manta')."/configManta.py", $pars, true);
$parser->exec("python {$manta_folder}/runWorkflow.py", "--mode local --jobs {$threads} --memGb 4", false);

//copy files to output folder
$struc = $manta_folder."/results/variants/diploidSV.vcf.gz";
if (!empty($t_bam)) $struc = $manta_folder."/results/variants/somaticSV.vcf.gz";
$parser->exec("cp", "{$struc} {$out}", false);
$parser->exec("cp", "{$struc}.tbi {$out}.tbi", false);

$small = $manta_folder."/results/variants/candidateSmallIndels.vcf.gz";
if (!empty($smallIndels))
{
	$parser->exec("cp", "{$small} {$smallIndels}", false);
	$parser->exec("cp", "{$small}.tbi {$smallIndels}.tbi", false);
}