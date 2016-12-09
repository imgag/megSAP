<?php

/**
	@page somatic
	@todo run DNA and RNA in parallel
	@todo add abra flag, add indel realignment for all files
	@todo check if all files are available (RNA and DNA, e.g. fastq)
	@todo consider haplotyping for somatic variants and germline variants close to each other
	@todo add germline variants for ADME genes
*/

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";
require_once($basedir."Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_emed", "\$Rev: 922 $", "Differential analysis of tumor/reference exomes. Can combine DNA and RNA data and results are suitable for current vaccination projects.");
$parser->addString("p_folder", "Folder that contains Sample folders.", false);
$parser->addString("t_dna_id",  "Tumor sample BAM file.", false);
$parser->addString("n_dna_id",  "Reference sample BAM file.", false);
$parser->addString("t_rna_id",  "Tumor sample BAM file.", true, "na");
$parser->addString("t_rna_fo",  "Tumor sample BAM file.", true, "na");
$parser->addString("n_rna_id",  "Normal sample BAM file.", true, "na");
$parser->addString("n_rna_fo",  "Normal sample BAM file.", true, "na");
$parser->addString("o_folder", "Output folder.", false);
// optional
$parser->addInfile("sys_tum_dna",  "Tumor DNA processing system INI file (used during annotation, determined from 't_dna_id' by default).", true);
$parser->addInfile("sys_nor_dna",  "Normal DNA processing system INI file (determined from 'n_dna_id' by default).", true);
$parser->addInfile("sys_tum_rna",  "Tumor RNA processing system INI file (used during annotation, determined from 't_rna_id' by default).", true);
$parser->addInfile("sys_nor_rna",  "Normal RNA processing system INI file (determined from 'n_rna_id' by default).", true);
$steps_all = array("ma", "vc", "fu", "an", "db", "co");
$parser->addString("steps", "Comma-separated list of processing steps to perform.", true, implode(",", $steps_all));
$parser->addFlag("rna_only", "Process (fastq,ma,vc) only RNA.", true);
$parser->addFlag("dna_only", "Process (ma,vc,db) only DNA.", true);
$parser->addFlag("no_germline", "Do not add analyze germline (mainly ADME variant list).", true);
$parser->addFlag("freebayes", "Use freebayes for somatic variant calling (default for tumor only). Can be reasonable e.g. for normal sample that is contaiminated with tumor.");
//$parser->addFlag("nsc", "Skip sample correlation check.");
extract($parser->parse($argv));

// determine steps to perform
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

// get processing systems
$t_dna_sys = load_system($sys_tum_dna, $t_dna_id);
$n_dna_sys = load_system($sys_nor_dna, $n_dna_id);
if($t_dna_sys["name_short"] != $n_dna_sys["name_short"]) trigger_error ("System tumor '".$t_dna_sys["name_short"]."' and system reference '".$n_dna_sys["name_short"]."' are different!", E_USER_WARNING);
if($t_dna_sys["build"] != $n_dna_sys["build"]) trigger_error ("System build tumor '".$t_dna_sys["build"]."' and system reference '".$n_dna_sys["build"]."' are different!", E_USER_ERROR);

// files and folders
$p_folder = rtrim($p_folder,'/')."/";
$o_folder = rtrim($o_folder,'/')."/";
if(!is_dir($o_folder))	mkdir($o_folder, 0775, true);
$t_dna_bam = $p_folder."Sample_".$t_dna_id."/".$t_dna_id.".bam";
$n_dna_bam = $p_folder."Sample_".$n_dna_id."/".$n_dna_id.".bam";
$s_dna_ann = $o_folder.$t_dna_id."-".$n_dna_id.".GSvar";
$s_dna_vcf = $o_folder.$t_dna_id."-".$n_dna_id."_var_annotated.vcf.gz";
$t_rna_bam = $p_folder."Sample_".$t_rna_id."/".$t_rna_id.".bam";
$n_rna_bam = $p_folder."Sample_".$n_rna_id."/".$n_rna_id.".bam";

// (1) run somatic_dna
$available_steps = array("ma", "vc", "an", "db");	//steps that can be used for somatic_dna script
if(count($tmp_steps=array_intersect($available_steps,$steps))>0 && !$rna_only)
{
	// somatic_dna pipelne starts with map and ends with db
	$args = "-sys_tum $sys_tum_dna -sys_nor $sys_nor_dna -steps ".implode(",",$tmp_steps)." ";
	$args .= "-filter_set somatic -min_af 0.05 -reduce_variants_strelka --log ".$o_folder."somatic_emed_dna_".date('YmdHis',mktime()).".log ";
	if($freebayes)	$args .= "-freebayes ";
	$parser->execTool("php $basedir/Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $t_dna_id -n_id $n_dna_id -o_folder $o_folder $args");
}

// (2) run somatic_rna, consider single sample mode, important for next step (s. (3))
$rna = 0;	//0 = no RNA, 1 = tumor sample RNA, 2 = tumor and reference sample RNA
if($t_rna_id!="na" && $t_rna_fo!="na" && ($n_rna_id=="na" || $n_rna_fo=="na"))	$rna = 1;
if($t_rna_id!="na" && $t_rna_fo!="na" && $n_rna_id!="na" && $n_rna_fo!="na")	$rna = 2;
$available_steps = array("ma","fu","db");	//steps that can be used for somatic_rna script
if(count($tmp_steps=array_intersect($available_steps,$steps))>0 && !$dna_only)
{	
	if($rna>0)
	{
		if(!is_dir($o_folder))	mkdir($o_folder, 0775, true);

		// somatic_rna pipeline starts with map and ends with an
		$args = "";
		$t_rna_sys = load_system($sys_tum_rna, $t_rna_id);
		$args .= "-sys_tum $sys_tum_rna ";
		if($rna>1)
		{
			$n_rna_sys = load_system($sys_nor_rna, $n_rna_id);
			$args .= "-sys_nor $sys_nor_rna ";
		}
		$args .= "-t_folder $t_rna_fo -n_folder $n_rna_fo --log ".$o_folder."somatic_emed_rna_".date('YmdHis',mktime()).".log ";
		$args .= "-steps ".implode(",", $tmp_steps)." ";
		$parser->execTool("php $basedir/Pipelines/somatic_rna.php", "-p_folder $p_folder -t_id $t_rna_id -n_id $n_rna_id -o_folder $o_folder ".$args);

		// compare RNA tumor and DNA tumor (tumor and normal is compared by the pipelines itself
		$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in1 $t_dna_bam -in2 $t_rna_bam -bam -max_snps 4000", true);
		$parts = explode(":", $output[0][1]);
		if ($parts[1]<0.8)	trigger_error("The genotype correlation of samples $t_bam and $n_bam is ".$parts[1]."; it should be above 0.8!", E_USER_ERROR);
	}
}

// (3) combine results from DNA / RNA
$g_dna_vcf = $o_folder."/".$n_dna_id."_var_annotated_adme.vcf.gz";
if(in_array("co", $steps))
{
	// (3) annotate RNA 
	$vaf_options = " -depth"; #-ref ".get_path("local_data")."/hg19.fa";
	if ($rna!=0)
	{
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $s_dna_ann -bam $t_rna_bam -out $s_dna_ann -name rna_tum $vaf_options", true);
		if ($rna==2)
		{
			$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $s_dna_ann -bam $n_rna_bam -out $s_dna_ann -name rna_ref $vaf_options", true);
		}
	}

	// (3) determine tumor content
	$parser->exec(get_path("ngs-bits")."EstimateTumorContent", "-tu $s_dna_ann -tu_bam ".$t_dna_bam." -no_bam ".$n_dna_bam, true);
	
	if(!$no_germline)
	{
		// (4) germline variants ADME genes
		$target_adme = "/mnt/projects/research/eMed-HCC/+IKP/ADME_Genes_hg19_eMED_20161206.bed";
		$tmp_folder = $parser->tempFolder("analyze");
		$adme_vcffile = $tmp_folder."/".$n_dna_id.".vcf.gz";
		$adme_varfile = $tmp_folder."/".$n_dna_id.".GSvar";
		$nor_bam = $p_folder."/Sample_".$n_dna_id."/".$n_dna_id.".bam";
		$extras_vc = "-target $target_adme ";
		$parser->execTool("php ".$basedir."NGS/vc_freebayes.php", "-bam $nor_bam -out $adme_vcffile -build ". $n_dna_sys['build']." $extras_vc");
		$parser->execTool("php ".$basedir."Pipelines/annotate.php", "-out_name $n_dna_id -out_folder $tmp_folder -system $sys_tum_dna -vcf $adme_vcffile");
		copy($tmp_folder."/".$n_dna_id."_var_annotated.vcf.gz",$g_dna_vcf);
		copy($tmp_folder."/".$n_dna_id.".GSvar",$o_folder."/".$n_dna_id."_adme.GSvar");
	}
}

// (5) prepare IGV-session for all files
if(is_dir($p_folder) && is_dir($o_folder))
{
	$rel_path = relative_path($o_folder, $p_folder);
	$igv_session = array();
	$igv_session[] = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>";
	$igv_session[] = "<Session genome=\"hg19\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"chr2:1544321-1544606\" path=\".\" version=\"8\">";
	$igv_session[] = "    <Resources>";
	if(is_file($s_dna_vcf))	$igv_session[] = "        <Resource path=\"".$rel_path."/Somatic_".basename($s_dna_vcf,".vcf.gz")."/".basename($s_dna_vcf)."\"/>";
	if(is_file($s_dna_vcf))	$igv_session[] = "        <Resource path=\"".$rel_path."/Somatic_".basename($g_dna_vcf,".vcf.gz")."/".basename($g_dna_vcf)."\"/>";
	if(is_file($t_dna_bam))	$igv_session[] = "        <Resource path=\"".$rel_path."/Sample_".basename($t_dna_bam,".bam")."/".basename($t_dna_bam)."\"/>";
	if(is_file($n_dna_bam))	$igv_session[] = "        <Resource path=\"".$rel_path."/Sample_".basename($n_dna_bam,".bam")."/".basename($n_dna_bam)."\"/>";
	if(is_file($t_rna_bam))	$igv_session[] = "        <Resource path=\"".$rel_path."/Sample_".basename($t_rna_bam,".bam")."/".basename($t_rna_bam)."\"/>";
	if(is_file($n_rna_bam))	$igv_session[] = "        <Resource path=\"".$rel_path."/Sample_".basename($n_rna_bam,".bam")."/".basename($n_rna_bam)."\"/>";
	$igv_session[] = "    </Resources>";
	$igv_session[] = "    <HiddenAttributes>";
	$igv_session[] = "        <Attribute name=\"NAME\"/>";
	$igv_session[] = "        <Attribute name=\"DATA FILE\"/>";
	$igv_session[] = "        <Attribute name=\"DATA TYPE\"/>";
	$igv_session[] = "    </HiddenAttributes>";
	$igv_session[] = "</Session>";
	file_put_contents($o_folder."/".$t_dna_id."-".$n_dna_id."_igv.xml",$igv_session);
}
else
{
	trigger_error("IGV-Session-File was not created. Folder $p_folder or $o_folder does not exist.",E_USER_WARNING);
}