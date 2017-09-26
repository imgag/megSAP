<?php

/**
	@page somatic_emed

	@todo run DNA and RNA in parallel
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_emed", "Differential analysis of tumor/reference exomes. Can combine DNA and RNA data and results are suitable for current vaccination projects.");
$parser->addString("p_folder", "Folder that contains Sample folders.", false);
$parser->addString("t_dna_id",  "Tumor sample DNA processing ID.", false);
$parser->addString("n_dna_id",  "Reference sample DNA processing ID.", false);
$parser->addString("t_rna_id",  "Tumor sample RNA processing ID.", true, "na");
$parser->addString("t_rna_fo",  "Tumor sample RNA data folder.", true, "na");
$parser->addString("n_rna_id",  "Normal sample RNA processing ID.", true, "na");
$parser->addString("n_rna_fo",  "Normal sample RNA data folder.", true, "na");
$parser->addString("o_folder", "Output folder.", false);
// optional
$parser->addInfile("t_dna_sys",  "Tumor DNA processing system INI file (used during annotation, determined from 't_dna_id' by default).", true);
$parser->addInfile("n_dna_sys",  "Normal DNA processing system INI file (determined from 'n_dna_id' by default).", true);
$parser->addInfile("t_rna_sys",  "Tumor RNA processing system INI file (used during annotation, determined from 't_rna_id' by default).", true);
$parser->addInfile("n_rna_sys",  "Normal RNA processing system INI file (determined from 'n_rna_id' by default).", true);
$steps_all = array("ma", "vc", "fu", "an", "db", "co");
$parser->addString("steps", "Comma-separated list of processing steps to perform.", true, implode(",", $steps_all));
$parser->addFlag("rna_only", "Process (fastq,ma,vc) only RNA.", true);
$parser->addFlag("dna_only", "Process (ma,vc,db) only DNA.", true);
$parser->addFlag("no_germline", "Do not analyze germline variants ('ADME gene variant list').", true);
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
$t_dna_sys_file = $t_dna_sys;
$t_dna_sys = load_system($t_dna_sys_file, $t_dna_id);
$n_dna_sys_file = $n_dna_sys;
$n_dna_sys = load_system($n_dna_sys_file, $n_dna_id);
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
$s_dna_seg = $o_folder.$t_dna_id."-".$n_dna_id."_cnvs.seg";
$t_rna_bam = $p_folder."Sample_".$t_rna_id."/".$t_rna_id.".bam";
$n_rna_bam = $p_folder."Sample_".$n_rna_id."/".$n_rna_id.".bam";

// (1) run somatic_dna
$available_steps = array("ma", "vc", "an", "db");	//steps that can be used for somatic_dna script
if(count($tmp_steps=array_intersect($available_steps,$steps))>0 && !$rna_only)
{
	// somatic_dna pipelne starts with map and ends with db
	$args = "-t_sys $t_dna_sys_file -n_sys $n_dna_sys_file -steps ".implode(",",$tmp_steps)." ";
	$args .= "-filter_set not-coding-splicing ";
	$args .= "-min_af 0.05 --log ".$o_folder."somatic_emed_dna_".date('YmdHis',mktime()).".log ";
	if($freebayes)	$args .= "-freebayes ";
	$parser->execTool("Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $t_dna_id -n_id $n_dna_id -o_folder $o_folder $args");
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
		$t_rna_sys_file = $t_rna_sys;
		$t_rna_sys = load_system($t_rna_sys_file, $t_rna_id);
		$args .= "-t_sys $t_rna_sys_file ";
		if($rna>1)
		{
			$n_rna_sys_file = $n_rna_sys;
			$n_rna_sys = load_system($n_rna_sys_file, $n_rna_id);
			$args .= "-n_sys $n_rna_sys_file ";
		}
		$args .= "-t_folder $t_rna_fo -n_folder $n_rna_fo --log ".$o_folder."somatic_emed_rna_".date('YmdHis',mktime()).".log ";
		$args .= "-steps ".implode(",", $tmp_steps)." ";
		$parser->execTool("Pipelines/somatic_rna.php", "-p_folder $p_folder -t_id $t_rna_id -n_id $n_rna_id -o_folder $o_folder ".$args);

		// compare RNA tumor and DNA tumor (tumor and normal is compared by the pipelines itself
		$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in1 $t_dna_bam -in2 $t_rna_bam -bam -max_snps 4000", true);
		$parts = explode(":", $output[0][1]);
		if ($parts[1]<0.8)	trigger_error("The genotype correlation of samples $t_bam and $n_bam is ".$parts[1]."; it should be above 0.8!", E_USER_ERROR);
	}
}

// (3) combine results from DNA / RNA
$g_dna_vcf = $o_folder."/".$n_dna_id."_adme_var_annotated.vcf.gz";
if(in_array("co", $steps))
{
	// (3) annotate RNA 
	if ($rna!=0)
	{
		$vaf_options = " -depth";
		
		$tmp = $parser->tempFile(".vcf");
		$parser->exec("bgzip", "-dc $s_dna_vcf > $tmp", false);	// no output logging, because Toolbase::extractVersion() does not return
		
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $tmp -bam $t_rna_bam -out $tmp -name rna_tum $vaf_options", true);
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $s_dna_ann -bam $t_rna_bam -out $s_dna_ann -name rna_tum $vaf_options", true);
		if ($rna==2)
		{
			$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $tmp -bam $n_rna_bam -out $tmp -name rna_ref $vaf_options", true);
			$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $s_dna_ann -bam $n_rna_bam -out $s_dna_ann -name rna_ref $vaf_options", true);
		}
		
		$parser->exec("bgzip", "-cf $tmp > $s_dna_vcf", false);	// no output logging, because Toolbase::extractVersion() does not return
		$parser->exec("tabix", "-fp vcf $s_dna_vcf", false);	// no output logging, because Toolbase::extractVersion() does not return
	}

	// (4) germline variants ADME genes
	if(!$no_germline)
	{
		$target_adme = "/mnt/projects/research/eMed-HCC/+IKP/ADME_Genes_hg19_eMED_20161206.bed";
		$tmp_folder = $parser->tempFolder("analyze");
		$adme_vcffile = $tmp_folder."/".$n_dna_id.".vcf.gz";
		$adme_varfile = $tmp_folder."/".$n_dna_id.".GSvar";
		$nor_bam = $p_folder."/Sample_".$n_dna_id."/".$n_dna_id.".bam";
		$extras_vc = "-target $target_adme ";
		$parser->execTool("NGS/vc_freebayes.php", "-bam $nor_bam -out $adme_vcffile -build ". $n_dna_sys['build']." $extras_vc");
		$parser->execTool("Pipelines/annotate.php", "-out_name $n_dna_id -out_folder $tmp_folder -system $t_dna_sys_file -vcf $adme_vcffile");
		copy($tmp_folder."/".$n_dna_id."_var_annotated.vcf.gz",$g_dna_vcf);
		copy($tmp_folder."/".$n_dna_id.".GSvar",$o_folder."/".$n_dna_id."_adme.GSvar");
	}
}

// (5) prepare IGV-session for all files
if(is_dir($p_folder) && is_dir($o_folder))
{
	$rel_path = relative_path($o_folder, $p_folder);
	$igv_session = array();
	$igv_session[] = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
	$igv_session[] = "<Session genome=\"1kg_v37\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"all\" path=\".\" version=\"8\">\n";
	$igv_session[] = "\t<Resources>";
	if(is_file($s_dna_seg))   $igv_session[] = "\t\t<Resource path=\"".$rel_path."/".$o_folder."/".basename($s_dna_seg)."\"/>\n";
	if(is_file($s_dna_vcf))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/".$o_folder."/".basename($s_dna_vcf)."\"/>\n";
	if(is_file($g_dna_vcf))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/".$o_folder."/".basename($g_dna_vcf)."\"/>\n";
	if(is_file($t_dna_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($t_dna_bam,".bam")."/".basename($t_dna_bam)."\"/>\n";
	if(is_file($n_dna_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($n_dna_bam,".bam")."/".basename($n_dna_bam)."\"/>\n";
	if(is_file($t_rna_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($t_rna_bam,".bam")."/".basename($t_rna_bam)."\"/>\n";
	if(is_file($n_rna_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($n_rna_bam,".bam")."/".basename($n_rna_bam)."\"/>\n";
	$igv_session[] = "\t</Resources>\n";
	$igv_session[] = "\t<HiddenAttributes>\n";
	$igv_session[] = "\t\t<Attribute name=\"NAME\"/>\n";
	$igv_session[] = "\t\t<Attribute name=\"DATA FILE\"/>\n";
	$igv_session[] = "\t\t<Attribute name=\"DATA TYPE\"/>\n";
	$igv_session[] = "\t</HiddenAttributes>\n";
	$igv_session[] = "</Session>";
	file_put_contents($o_folder."/".$t_dna_id."-".$n_dna_id."_igv.xml",$igv_session);
}
else
{
	trigger_error("IGV-Session-File was not created. Folder $p_folder or $o_folder does not exist.",E_USER_WARNING);
}