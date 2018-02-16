<?php

/**
	@page somatic_dna_rna
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_dna_rna", "Differential analysis of tumor/reference exomes and transcriptomes.");
$parser->addString("p_folder", "Folder that contains sample folders.", false);

$parser->addString("t_dna_id",  "Tumor sample DNA processing ID.", false);
$parser->addString("n_dna_id",  "Reference sample DNA processing ID.", false);

$parser->addString("t_rna_id",  "Tumor sample RNA processing ID.", true, "na");
$parser->addString("n_rna_id",  "Normal sample RNA processing ID.", true, "na");

// optional
$parser->addString("process", "Enable processing for DNA pair, RNA pair, germline analysis for DNA normal sample and combination steps.", true, "dna,rna,germline,co,igv");
$parser->addString("steps_dna", "Processing steps for DNA samples (somatic_dna).", true, "vc,an,db");
$parser->addString("steps_rna", "Processing steps for RNA samples (somatic_rna).", true, "rc,fu");

$parser->addString("germline_preset", "Germline analysis preset:\ndefault: use target germline_target (if specified) for analysis of normal DNA\nnearby: run nearby-germline-variants strategy", true, "default");
$parser->addString("germline_target", "Target region for default germline analysis.", true, "");
$parser->addString("germline_suffix", "File name suffix to append to germline analysis results.", true, "");

$parser->addString("o_folder", "Output folder, defaults to <p_folder>/Somatic_<t_dna_id>-<n_dna_id>.", true, "default");

$parser->addInfile("t_dna_sys",  "Tumor DNA processing system INI file (used during annotation, determined from 't_dna_id' by default).", true);
$parser->addInfile("n_dna_sys",  "Normal DNA processing system INI file (determined from 'n_dna_id' by default).", true);
$parser->addInfile("t_rna_sys",  "Tumor RNA processing system INI file (used during annotation, determined from 't_rna_id' by default).", true);
$parser->addInfile("n_rna_sys",  "Normal RNA processing system INI file (determined from 'n_rna_id' by default).", true);

$parser->addString("filter_set", "Filter set for somatic variants (somatic_dna).", true, "not-coding-splicing");
$parser->addFloat("min_af", "Allele frequency detection cut-off (somatic_dna).", true, 0.05);
// TODO add further somatic_dna options, e.g. freebayes (?)

$parser->addFlag("skip_correlation", "Skip sample correlation check.");

// extract command line arguments
extract($parser->parse($argv));
$process = explode(",", $process);
//$steps_dna = explode(",", $steps_dna);
//$steps_rna = explode(",", $steps_rna);
$rna_tum_available = $t_rna_id !== "na";
$rna_ref_available = $n_rna_id !== "na";

// log version
$parser->log("Pipeline revision: " . repository_revision(true));

// resolve default output folder
if ($o_folder === "default")
{
	$o_folder = "{$p_folder}/Somatic_{$t_dna_id}-{$n_dna_id}";
}

// create output folder if necessary
if (!file_exists($o_folder))
{
	mkdir($o_folder);
	if (!chmod($o_folder, 0777))
	{
		trigger_error("Could not change privileges of folder '{$o_folder}'!", E_USER_ERROR);
	}
}

// extract systems
$t_dna_system = load_system($t_dna_sys, $t_dna_id);
$n_dna_system = load_system($n_dna_sys, $n_dna_id);

$t_rna_system = $t_rna_id === "na" ? "na" : load_system($t_rna_sys, $t_rna_id);
$n_rna_system = $n_rna_id === "na" ? "na" : load_system($n_rna_sys, $n_rna_id);

// define file paths

// bam files
$t_dna_bam = "{$p_folder}/Sample_{$t_dna_id}/{$t_dna_id}.bam";
$n_dna_bam = "{$p_folder}/Sample_{$n_dna_id}/{$n_dna_id}.bam";
$t_rna_bam = "{$p_folder}/Sample_{$t_rna_id}/{$t_rna_id}.bam";
$n_rna_bam = "{$p_folder}/Sample_{$n_rna_id}/{$n_rna_id}.bam";

// somatic variants
$som_dna_vcf = "{$o_folder}/{$t_dna_id}-{$n_dna_id}_var_annotated.vcf.gz";
$som_dna_gsvar = "{$o_folder}/{$t_dna_id}-{$n_dna_id}.GSvar";
$som_dna_seg = "{$o_folder}/{$t_dna_id}-{$n_dna_id}_cnvs.seg";

// germline variants, placed in somatic folder
$germline_dna_vcf = "{$o_folder}/{$n_dna_id}{$germline_suffix}_var_annotated.vcf.gz";
$germline_dna_gsvar = "{$o_folder}/{$n_dna_id}{$germline_suffix}.GSvar";

// run somatic_dna
if (in_array("dna", $process))
{
	$somatic_dna_args = [
		"-p_folder", $p_folder,
		"-t_id", $t_dna_id,
		"-n_id", $n_dna_id,
		"-o_folder", $o_folder,
		"-steps", $steps_dna,
		"-filter_set", $filter_set,
		"-min_af", $min_af,
		"--log", "{$o_folder}/somatic_dna_" . date('YmdHis', mktime()) . ".log"
	];
	if (isset($t_dna_sys))
	{
		$somatic_dna_args[] = "-t_sys {$t_dna_sys}";
	}
	if (isset($n_dna_sys))
	{
		$somatic_dna_args[] = "-n_sys {$n_dna_sys}";
	}
	if ($skip_correlation)
	{
		$somatic_dna_args[] = "-nsc";
	}

	$parser->execTool("Pipelines/somatic_dna.php", implode(" ", $somatic_dna_args));
}

// run somatic_rna
if (in_array("rna", $process) && $rna_tum_available)
{
	$somatic_rna_args = [
		"-p_folder", $p_folder,
		"-t_id", $t_rna_id,
		"-n_id", $n_rna_id,
		"-o_folder", $o_folder,
		"-steps", $steps_rna,
		"--log", "{$o_folder}/somatic_rna_" . date('YmdHis', mktime()) . ".log"
		];
	if (isset($t_rna_sys))
	{
		$somatic_rna_args[] = "-t_sys {$t_rna_sys}";
	}
	if (isset($n_rna_sys))
	{
		$somatic_rna_args[] = "-n_sys {$n_rna_sys}";
	}
	if ($skip_correlation)
	{
		$somatic_rna_args[] = "-skip_correlation";
	}
	$parser->execTool("Pipelines/somatic_rna.php", implode(" ", $somatic_rna_args));
}

// run germline analysis on normal DNA sample
if (in_array("germline", $process))
{
	if ($germline_preset === "default")
	{
		$germline_dna_tmp = $parser->tempFolder("germline_dna");
		$germline_dna_vcf_tmp = "{$germline_dna_tmp}/{$n_dna_id}_var_annotated.vcf.gz";
		$germline_dna_gsvar_tmp = "{$germline_dna_tmp}/{$n_dna_id}.GSvar";
		
		$vc_args = [
			"-bam", $n_dna_bam,
			"-out", $germline_dna_vcf_tmp,
			"-build", $n_dna_system['build']
		];
		
		if ($germline_target !== "")
		{
			$vc_args[] = "-target {$germline_target}";
		}
		
		$parser->execTool("NGS/vc_freebayes.php", implode(" ", $vc_args));
		$parser->execTool("Pipelines/annotate.php", "-out_name {$n_dna_id} -out_folder {$germline_dna_tmp} -system {$n_dna_sys} -vcf {$germline_dna_vcf_tmp}");
		
		$parser->copyFile($germline_dna_vcf_tmp, $germline_dna_vcf);
		$parser->copyFile($germline_dna_gsvar_tmp, $germline_dna_gsvar);
	}
	else if ($germline_preset === "nearby")
	{
		$parser->execTool("NGS/somatic_nearby_germline_variants.php", "-somatic_var {$som_dna_gsvar} -normal_bam {$n_dna_bam} -n_dna_id {$n_dna_id} -n_dna_sys {$n_dna_sys}");
	}
}

// combine somatic DNA variants (vcf, GSvar) with RNA frequency/depth
if (in_array("co", $process) && ($rna_tum_available || $rna_ref_available))
{
	// check correlation of RNA tumor and DNA tumor, warn if disabled
	if ($skip_correlation && $rna_tum_available)
	{
		trigger_error("Sample correlation check tumor DNA / tumor RNA has been disabled!", E_USER_WARNING);
	}
	else if ($rna_tum_available)
	{
		$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in {$t_dna_bam} {$t_rna_bam} -mode bam -max_snps 4000", true);
		$correlation = explode("\t", $output[0][1])[3];
		if ($correlation < 0.8)
		{
			trigger_error("The genotype correlation of samples {$t_dna_bam} and {$t_rna_bam} is {$correlation}, it should be at least 0.8!", E_USER_ERROR);
		}
	}

	// annotate somatic variants (from exomes) with RNA frequency/depth

	// decompress somatic vcf
	$som_dna_vcf_tmp = $parser->tempFile("var_annotated.vcf");
	$parser->exec("bgzip", "-dc {$som_dna_vcf} > {$som_dna_vcf_tmp}", false);

	if ($rna_tum_available)
	{
		// annotate vcf
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency",
			"-in {$som_dna_vcf_tmp} -bam {$t_rna_bam} -out {$som_dna_vcf_tmp} -name rna_tum -depth", true);
		// annotate GSvar
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency",
			"-in {$som_dna_gsvar} -bam {$t_rna_bam} -out {$som_dna_gsvar} -name rna_tum -depth", true);
	}

	if ($rna_ref_available)
	{
		// annotate vcf
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency",
			"-in {$som_dna_vcf_tmp} -bam {$n_rna_bam} -out {$som_dna_vcf_tmp} -name rna_ref -depth", true);
		// annotate GSvar
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency",
			"-in {$som_dna_gsvar} -bam {$n_rna_bam} -out {$som_dna_gsvar} -name rna_ref -depth", true);
	}

	// compress + index somatic vcf
	$parser->exec("bgzip", "-cf {$som_dna_vcf_tmp} > {$som_dna_vcf}", false);
	$parser->exec("tabix", "-fp vcf {$som_dna_vcf}", false);
}

// create IGV session file
if (in_array("igv", $process))
{
	$igv_session_file = "{$o_folder}/{$t_dna_id}-{$n_dna_id}_igv.xml";
	$igv_tracks = implode(" ", array_filter([
		$som_dna_seg,
		$som_dna_vcf,
		$germline_dna_vcf,
		$t_dna_bam,
		$n_dna_bam,
		$t_rna_bam,
		$n_rna_bam
		], "file_exists"));
	$parser->execTool("NGS/igv_session.php", "-out {$igv_session_file} -in {$igv_tracks} -relative");
}