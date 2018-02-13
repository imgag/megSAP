<?php

/*
	@page somatic_rna
 */

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";
require_once($basedir."Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// command line arguments
$parser = new ToolBase("somatic_rna", "Differential analysis of tumor-normal RNA sample pairs.");
$parser->addString("p_folder","Folder containing sample subfolders with fastqs (Sample_GSXYZ).",false);
$parser->addString("t_id", "Tumor sample processing-ID (e.g. GSxyz_01). There should be a folder 'Sample_tsid' that contains fastq-files marked with the id within the p_folder.", false);

// optional arguments
$parser->addString("n_id", "Normal sample processing-ID (e.g. GSxyz_01). There should be a folder 'Sample_nsid' that contains fastq-files marked with the id within the p_folder.", true, "na");
$parser->addString("o_folder", "Output folder, defaults to <p_folder>/Somatic_<t_id>-<n_id>.", true, "default");

$steps_all = array("ma", "rc", "fu", "vc", "an");
$parser->addString("steps", "Comma-separated list of processing steps to perform:\n" .
	"ma = mapping, rc = compare read counts, fu = find somatic fusions, vc = variant calling (strelka), an = annotation of variants", true, "rc,fu,vc,an");

$parser->addFlag("skip_correlation", "Skip sample correlation check.");

$parser->addInfile("t_sys",  "Tumor processing system INI file (determined from 't_id' by default).", true);
$parser->addInfile("n_sys",  "Reference processing system INI file (determined from 'n_id' by default).", true);

$parser->addFlag("debug", "Enable debug mode for testing purposes.");

// extract command line arguments
extract($parser->parse($argv));
$steps = explode(",", $steps);

// log version
$parser->log("Pipeline revision: " . repository_revision(true));

// extract systems
$t_system = load_system($t_sys, $t_id);
$n_system = $n_id === "na" ? "na" : load_system($n_sys, $n_id);

// if only tumor specified, completely hand off to analyze_rna.php
if ($n_id == "na")
{
	if (in_array("ma", $steps))
	{
		$parser->execTool("Pipelines/analyze_rna.php",
			"-folder {$p_folder}/Sample_{$t_id} -name {$t_id} -system {$t_sys} -steps ma,rc,an,fu,db --log {$p_folder}/Sample_{$t_id}/analyze_".date('YmdHis',mktime()).".log");
	}
	exit(0);
}

// resolve default output folder
if ($o_folder === "default")
{
	$o_folder = "{$p_folder}/Somatic_{$t_id}-{$n_id}";
}

// warn about different system
if ($t_system["name_short"] !== $n_system["name_short"])
{
	trigger_error("Tumor and normal sample have different processing systems!", E_USER_WARNING);
}

// fail on different references
if ($t_system["build"] !== $n_system["build"])
{
	trigger_error("Tumor and normal sample have different genome builds!", E_USER_ERROR);
}

// mapping
if (in_array("ma", $steps))
{
	$parser->execTool("Pipelines/analyze_rna.php",
		"-folder Sample_{$t_id} -name {$t_id} -system {$t_sys} -steps ma,rc,an,fu,db --log Sample_{$t_id}/analyze_".date('YmdHis',mktime()).".log");
	$parser->execTool("Pipelines/analyze_rna.php",
		"-folder Sample_{$n_id} -name {$n_id} -system {$n_sys} -steps ma,rc,an,fu,db --log Sample_{$n_id}/analyze_".date('YmdHis',mktime()).".log");
}

// define file paths

// bam files
$tum_bam = "{$p_folder}/Sample_{$t_id}/{$t_id}.bam";
$ref_bam = "{$p_folder}/Sample_{$n_id}/{$n_id}.bam";

// read count files
$tum_counts = "{$p_folder}/Sample_{$t_id}/{$t_id}_counts.tsv";
$ref_counts = "{$p_folder}/Sample_{$n_id}/{$n_id}_counts.tsv";
$som_counts = "{$o_folder}/{$t_id}-{$n_id}_counts.tsv";

// fusion files
$tum_fusions = "{$p_folder}/Sample_{$t_id}/{$t_id}_var_fusions.tsv";
$ref_fusions = "{$p_folder}/Sample_{$n_id}/{$n_id}_var_fusions.tsv";
$som_fusions = "{$o_folder}/{$t_id}-{$n_id}_var_fusions.tsv";

// strelka files
$som_variants = "{$o_folder}/{$t_id}-{$n_id}_var.vcf.gz";
$som_variants_annotated = "{$o_folder}/{$t_id}-{$n_id}_var_annotated.vcf.gz";

// annotation
$som_gsvar = "{$o_folder}/{$t_id}-{$n_id}.GSvar";

// check sample correlation, warn if disabled
if ($skip_correlation)
{
	trigger_error("Sample correlation check has been disabled!", E_USER_WARNING);
}
else
{
	$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in {$tum_bam} {$ref_bam} -mode bam -max_snps 4000", true);
	$correlation = explode("\t", $output[0][1])[3];
	if ($correlation < 0.8)
	{
		trigger_error("The genotype correlation of samples {$tum_bam} and {$ref_bam} is {$correlation}, it should be at least 0.8!", E_USER_ERROR);
	}
}

// compare read counts
if (in_array("rc", $steps))
{
	trigger_error("Running step 'rc' ...", E_USER_NOTICE);
	$parser->execTool("NGS/rc_compare.php", "-in1 {$tum_counts} -in2 {$ref_counts} -out {$som_counts}");
	$parser->execTool("NGS/rc_annotate.php", "-in {$som_counts} -out {$som_counts}");
}

// find tumor sample fusions, which are not in normal sample fusions
if (in_array("fu", $steps))
{
	trigger_error("Running step 'fu' ...", E_USER_NOTICE);
	
	$fusions1 = Matrix::fromTSV($tum_fusions);	//tumor
	$idx_tleft = $fusions1->getColumnIndex("LeftBreakpoint");
	$idx_tright = $fusions1->getColumnIndex("RightBreakpoint");
	$fusions2 = Matrix::fromTSV($ref_fusions);	//normal
	$idx_nleft = $fusions2->getColumnIndex("LeftBreakpoint");
	$idx_nright = $fusions2->getColumnIndex("RightBreakpoint");

	$fusions_somatic = new Matrix();
	$fusions_somatic->setHeaders($fusions1->getHeaders());
	for($i=0;$i<$fusions1->rows();++$i)
	{
		$somatic = true;

		$r_tum = $fusions1->getRow($i);

		for($j=0;$j<$fusions2->rows();++$j)
		{
			$r_nor = $fusions2->getRow($j);

			if($r_tum[$idx_tleft]==$r_nor[$idx_nleft] && $r_tum[$idx_tright]==$r_nor[$idx_nright])	$somatic = false;
		}

		if($somatic)	$fusions_somatic->addRow($r_tum);
	}
	$fusions_somatic->toTSV($som_fusions);
}

// run strelka
if (in_array("vc", $steps))
{
	trigger_error("Running step 'vc' ...", E_USER_NOTICE);
	if ($debug)
	{
		$debug_args = " -debug_region chr22";
	}
	$parser->execTool("NGS/vc_strelka2.php", "-t_bam {$tum_bam} -n_bam {$ref_bam} -out {$som_variants}{$debug_args}");
}

// annotate strelka results
if (in_array("an", $steps))
{
	trigger_error("Running step 'an' ...", E_USER_NOTICE);
	
	// annotate vcf
	$parser->execTool("Pipelines/annotate.php", "-out_name {$t_id}-{$n_id} -out_folder {$o_folder} -system {$t_sys} -vcf {$som_variants} -t_col {$t_id} -n_col {$n_id}");
	
	// convert vcf to GSvar
	$parser->execTool("NGS/vcf2gsvar_somatic.php", "-in {$som_variants_annotated} -out {$som_gsvar} -t_col {$t_id}");
	
	// annotate GSvar (dbNFSP, NGSD, frequency + depth)
	$parser->execTool("NGS/an_dbNFSPgene.php", "-in {$som_gsvar} -out {$som_gsvar}");
	$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in {$som_gsvar} -out {$som_gsvar} -psname {$t_id} -mode somatic", true);
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in {$som_gsvar} -bam {$tum_bam} -out {$som_gsvar} -name rna_tum -depth", true);
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in {$som_gsvar} -bam {$ref_bam} -out {$som_gsvar} -name rna_ref -depth", true);
}