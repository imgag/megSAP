<?php

/**
  @page vc_manta
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_manta", "Call of germline/somatic structural variants with manta. Creates an VCF file.");
$parser->addOutfile("out", "Output VCF file (gzipped and tabix indexed).", false);
//optional
$parser->addInfileArray("bam", "Normal BAM file(s). Only one normal BAM file allowed for somatic mode.", true, false);
$parser->addInfile("t_bam", "Tumor BAM file, for somatic mode.", true, false);
$parser->addOutfile("smallIndels", "Output VCF file for candidate small indels (gzipped and tabix indexed).", true);
$parser->addString("evid_dir", "Output folder for BAM files containing evidence reads.",true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addInfile("target",  "Enrichment target BED file (used for flagging off-target variants).", true);
$parser->addFlag("exome", "If set, manta settings for exome/panel analysis are used (no depth filtering).");
$parser->addFlag("rna", "If set, manta settings for RNA fusion calling are used.");
$parser->addInt("threads", "Number of threads used.", true, 4);
//debugging options
$parser->addEnum("config_preset", "Use preset configuration.", true, array("default", "high_sensitivity"), "default");
$parser->addStringArray("regions", "Limit analysis to specified regions (for debugging).", true);
$parser->addString("temp", "Temporary folder for manta analysis (for debugging).", true, "auto");
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
$mode_germline = !isset($t_bam) && isset($bam) && !$rna;

// update environment
$python_bin = get_path("python27");
putenv("PATH=".dirname($python_bin).":".getenv("PATH"));

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

$args = [
	"--referenceFasta ".genome_fasta($build),
	"--runDir ".$manta_folder,
	"--config ".$config,
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
if ($mode_somatic || $mode_germline || $rna)
{
	array_push($args, "--normalBam", implode(" --normalBam ", $bam));
}
if (isset($regions))
{
	array_push($args, "--region", implode(" --region ", $regions));
}
if ($rna)
{
	array_push($args, "--rna");
}

//run manta
$parser->exec("{$python_bin} ".get_path('manta')."/configManta.py", implode(" ", $args), true);
$parser->exec("{$python_bin} {$manta_folder}/runWorkflow.py", "--mode local --jobs {$threads} --memGb ".(2*$threads), false);

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
else if ($rna)
{
	$outname = "rna";
}
$sv = "{$manta_folder}/results/variants/{$outname}SV.vcf.gz";

//combine BND of INVs to one INV in VCF
$sv_inv = "{$manta_folder}/results/variants/{$outname}SV_inv.vcf";
$parser->exec("{$python_bin} ".get_path('manta')."/../libexec/convertInversion.py", get_path("samtools")." ".genome_fasta($build)." {$sv} > {$sv_inv}");

//remove VCF lines with empty "REF". They are sometimes created from convertInversion.py but are not valid
$vcf_fixed = "{$temp_folder}/{$outname}SV_fixed.vcf";
$h = fopen2($sv_inv, "r");
$h2 = fopen2($vcf_fixed, "w");
while(!feof($h))
{
	$line = fgets($h);
	$parts = explode("\t", $line);
	if (count($parts)>3 && $parts[3]=="") continue;
	
	fputs($h2, $line);
}
fclose($h);
fclose($h2);

//sort variants
$vcf_sorted = "{$temp_folder}/{$outname}SV_sorted.vcf";
$parser->exec(get_path("ngs-bits")."VcfSort","-in {$vcf_fixed} -out $vcf_sorted", true);

// flag off-target variants
if (isset($target))
{
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

//Copy evidence bams in case of somatic/ tumor only sample
if(isset($evid_dir))
{
	create_directory($evid_dir);
	$parser->exec("cp","{$temp_folder}/mantaAnalysis/results/evidence/* {$evid_dir}",true);
}

$small = "{$manta_folder}/results/variants/candidateSmallIndels.vcf.gz";
if (isset($smallIndels))
{
	$parser->moveFile($small, $smallIndels);
	$parser->moveFile("{$small}.tbi", "{$smallIndels}.tbi");
}

?>