<?php

/**
  @page vc_manta
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_manta", "Call somatic structural variants with manta. Creates an VCF file.");
$parser->addOutfile("out", "Output file (gzipped and tabix indexed).", false);
//optional
$parser->addInfile("t_bam", "Tumor BAM file, for somatic mode.", true, false);
$parser->addInfileArray("bam", "Normal BAM file(s). Only one normal BAM file allowed for somatic mode.", true, false);
$parser->addFlag("exome", "If set, settings for exome analysis are used.", true);
$parser->addOutfile("smallIndels", "Output file for candidate small indels.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addString("temp", "Temporary folder for manta analysis.", true, "auto");
$parser->addStringArray("regions", "Limit analysis to specified regions.", true);
$parser->addInfile("target",  "Enrichment targets BED file.", true);

$parser->addEnum("config_preset", "Use preset configuration.", true, array("default", "high_sensitivity"), "default");

$parser->addInt("threads", "Number of threads used.", true, 4);

extract($parser->parse($argv));

// determine mode (somatic, tumor-only, germline)
if (!isset($t_bam) && !isset($bam))
{
	trigger_error("No input BAM file(s) specified!", E_USER_ERROR);
}
else if (isset($t_bam) && isset($bam) && count($bam) > 1)
{
	trigger_error("More than one normal sample specified for somatic analysis!", E_USER_ERROR);
}

$mode_somatic = isset($t_bam) && isset($bam) && count($bam) == 1;
$mode_tumor_only = isset($t_bam) && !isset($bam);
$mode_germline = !isset($t_bam) && isset($bam);

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

$args = [
	"--referenceFasta", $genome,
	"--runDir", $manta_folder,
	"--config", $config,
	"--outputContig",
	"--generateEvidenceBam"
];
if ($exome)
{
	array_push($args, "--exome");
}
if ($mode_somatic || $mode_tumor_only)
{
	array_push($args, "--tumorBam", $t_bam);
}
if ($mode_somatic || $mode_germline)
{
	array_push($args, "--normalBam", implode(" --normalBam ", $bam));
}
if (isset($regions)) {
	array_push($args, "--region", implode(" --region ", $regions));
}

//run manta
$parser->exec("python ".get_path('manta')."/configManta.py", implode(" ", $args), true);
$parser->exec("python {$manta_folder}/runWorkflow.py", "--mode local --jobs {$threads} --memGb 4", false);

//copy files to output folder
if ($mode_somatic)
{
	$outname = "somatic";
}
else if ($mode_tumor_only)
{
	$outname = "tumor";
}
else if ($mode_germline)
{
	$outname = "diploid";
}
$sv = "{$manta_folder}/results/variants/{$outname}SV.vcf.gz";

//sort variants
$vcf_sorted = "{$temp_folder}/{$outname}SV_sorted.vcf";
$parser->exec(get_path("ngs-bits")."VcfSort","-in $sv -out $vcf_sorted", true);

// flag off-target variants
if (isset($target))
{
	// TODO pre-process target?

	$vcf_filtered = "{$temp_folder}/{$outname}SV_filtered.vcf";
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $vcf_sorted -mark off-target -reg $target -out $vcf_filtered", true);
}
else
{
	$vcf_filtered = $vcf_sorted;
}

//zip and index output file
$parser->exec("bgzip", "-c $vcf_filtered > $out", true);
$parser->exec("tabix", "-p vcf $out", true);

$small = "{$manta_folder}/results/variants/candidateSmallIndels.vcf.gz";
if (isset($smallIndels))
{
	$parser->moveFile($small, $smallIndels);
	$parser->moveFile("{$small}.tbi", "{$smallIndels}.tbi");
}

?>