<?php

/**
	@page somatic
	@todo clip-overlap as default, CAVE: amplicon mode
	@todo cleanup fastq files
	@todo check if all files are available (fastq)
	@todo implement different steps
	@todo add indel realignment for all files
	@todo think about filter on-target (merging all target regions) and sample columns
*/

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";

require_once($basedir."Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("somatic_triplett", "\$Rev: 831 $", "Differential analysis of tumor/relapse/reference sample. DNA only.");
$parser->addString("p_folder","Folder containing sample subfolders with fastqs (Sample_GSXYZ).",false);
$parser->addString("t_dna_id",  "Tumor DNA processing ID.", false);
$parser->addString("r_dna_id",  "Relapse DNA processing ID.", false);
$parser->addString("n_dna_id",  "Normal DNA processing ID.", false);
$parser->addString("o_folder", "Output folder.", false);
//optional
$parser->addInfile("sys_tum", "Tumor processing system INI file (determined from 't_dna_id' by default).", true);
$parser->addInfile("sys_rel", "Relapse processing system INI file (determined from 'r_dna_id' by default).", true);
$parser->addInfile("sys_nor", "Normal processing system INI file (determined from 'n_dna_id' by default).", true);
$parser->addFlag("abra", "Use Abra for indel realignment.");
$steps_all = array("ma", "vc", "an", "combine");
$parser->addString("steps", "Comma-separated list of processing steps to perform.", true, implode(",", $steps_all));
$parser->addFlag("no_db_import", "Skip somatic variant import to db.", false);
$parser->addFlag("nsc", "Skip sample correlation check.");
extract($parser->parse($argv));

// determine steps to perform
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}
$t_dna_fo = $p_folder."/Sample_".$t_dna_id."/";
$r_dna_fo = $p_folder."/Sample_".$r_dna_id."/";
$n_dna_fo = $p_folder."/Sample_".$n_dna_id."/";
$o_folder .= "/";
if(!is_dir($t_dna_fo))	trigger_error("Could not find tumor folder '$t_dna_fo'.", E_USER_ERROR);
if(!is_dir($r_dna_fo))	trigger_error("Could not find relapse folder '$r_dna_fo'.", E_USER_ERROR);
if(!is_dir($n_dna_fo))	trigger_error("Could not find normal folder '$n_dna_fo'.", E_USER_ERROR);
if(!is_dir($o_folder))	mkdir($o_folder, 0775, true);

// run somatic_dna
//  map normal
if(in_array("ma", $steps))
{	
	$args = "-steps ma ";
	if(isset($sys_nor))	$args .= "-system $sys_nor ";
	$parser->execTool("php $basedir/Pipelines/analyze.php", "-folder ".$n_dna_fo." -name ".$n_dna_id." ".$args." --log ".$n_dna_fo."analyze_".date('YmdHis',mktime()).".log");
}

if(in_array("ma", $steps) || in_array("vc", $steps) || in_array("an", $steps))
{	
	// make somatic calling
	$tmp_steps = array();
	foreach($steps as $step)
	{
		if($step=="combine")	continue;
		$tmp_steps[] = $step;
	}
	if(!$no_db_import)	$tmp_steps[] = "db";
	$args = "-steps ".implode(",",$tmp_steps);
	if($nsc)	$args .= " -nsc";
	if($abra)	$args .= " -abra";

	//
	$tmp_sys = "";
	if(isset($sys_nor))	$tmp_sys .= "-sys_nor $sys_nor ";
	if(isset($sys_tum))	$tmp_sys .= "-sys_tum $sys_tum ";
	$parser->execTool("php $basedir/Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $t_dna_id -n_id $n_dna_id $tmp_sys -o_folder $o_folder $args -smn --log ".$o_folder."somatic_dna1_".date('YmdHis',mktime()).".log");

	//
	$tmp_sys = "";
	if(isset($sys_nor))	$tmp_sys .= "-sys_nor $sys_nor ";
	if(isset($sys_rel))	$tmp_sys .= "-sys_tum $sys_rel ";
	$parser->execTool("php $basedir/Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $r_dna_id -n_id $n_dna_id $tmp_sys -o_folder $o_folder $args -smn --log ".$o_folder."somatic_dna2_".date('YmdHis',mktime()).".log");
}

// combine results
if(in_array("combine", $steps))
{	
	$t_dna_bam = $t_dna_fo.$t_dna_id.".bam";
	$r_dna_bam = $r_dna_fo.$r_dna_id.".bam";
	$n_dna_bam = $n_dna_fo.$n_dna_id.".bam";
	$annotated1 = $o_folder.$t_dna_id."-".$n_dna_id.".GSvar";
	$annotated2 = $o_folder.$r_dna_id."-".$n_dna_id.".GSvar";
	$overview = $o_folder.$t_dna_id."-".$r_dna_id."-".$n_dna_id."_overview.tsv";
	$parser->exec(get_path("ngs-bits")."SampleOverview", "-in $annotated1 $annotated2 -add_cols interpro,Pathway_KEGG_full,Function_description,Pathway_BioCarta_full,ExAC,COSMIC -out ".$overview, true);
	$vaf_options = " -depth";
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $overview -bam $t_dna_bam -out $overview -name tum_".$t_dna_id." $vaf_options", true);
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $overview -bam $r_dna_bam -out $overview -name rel_".$r_dna_id." $vaf_options", true);
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $overview -bam $n_dna_bam -out $overview -name nor_".$n_dna_id." $vaf_options", true);
	$sys = load_system($sys_tum, $t_dna_id);
	$parser->exec("php $basedir/NGS/filter_tsv.php", "-in $overview -out $overview -type coding,non_synonymous -roi ".$sys['target_file'], true);
}

?>