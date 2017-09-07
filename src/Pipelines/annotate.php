<?php

/**
	@page annotate
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("annotate", "Annotate variants.");
$parser->addString("out_name", "Processed sample ID (e.g. 'GS120001_01').", false);
$parser->addString("out_folder", "Output folder.", false);
//optional
$parser->addInfile("system", "Processing system INI file (determined from 'out_name' by default).", true);
$parser->addString("t_col", "Column name of tumor sample (for Strelka VCF files only).", true, "na");
$parser->addString("n_col", "Column name of normal sample (for Strelka VCF files only).", true, "na");
$parser->addString("vcf", "Path to (bgzipped) VCF file (if different from {output_folder}/{out_name}_var.vcf.gz).", true, "");
$parser->addInt("thres", "Number of flanking bases around exons which are considered as splicing region.", true, 20);
$parser->addFlag("no_fc", "No format check (vcf/tsv).");
$parser->addFlag("no_ngsd", "No annotation of variants with information from NGSD.");
$parser->addFlag("multi", "Enable multi-sample mode.");
extract($parser->parse($argv));

//input file names
$vcffile = $out_folder."/".$out_name."_var.vcf.gz";
if(!empty($vcf))
{
	$vcffile = $vcf;
}

//unzip input VCF if necessary
if (ends_with($vcffile, ".vcf"))
{
	$vcf_unzipped = $vcffile;
}
else
{
	$vcf_unzipped = $parser->tempFile("_unzipped.vcf");
	$parser->exec("bgzip", "-cd $vcffile > $vcf_unzipped", false); //no output logging, because Toolbase::extractVersion() does not return
}

//output file names
$annfile = $parser->tempFile(".vcf");
$annfile_zipped = $out_folder."/".$out_name."_var_annotated.vcf.gz";	
$varfile = $out_folder."/".$out_name.".GSvar";	
$varfile_full = $out_folder."/".$out_name."_full.GSvar";
$stafile = $out_folder."/".$out_name."_stats_vc.qcML";

//get system
$sys = load_system($system, $out_name);
if ($sys['build']!="hg19" && $sys['build']!="GRCh37" && $sys['build']!="mm10")
{
	trigger_error("Unknown genome build ".$sys['build']." cannot be annotated!", E_USER_ERROR);
}

//annotate VCF
$args = array("-in $vcf_unzipped", "-thres $thres", "-out $annfile");
if ($sys['type']!="WGS")
{
	$args[] = "-no_updown";
}
$args[] = "-build ".$sys['build'];
$parser->execTool("NGS/an_snpeff.php", implode(" ", $args));

//check vcf file
if(!$no_fc)
{
	$parser->execTool("NGS/check_vcf.php", "-in $annfile");
}

//convert to GSvar file
if ($t_col=="na") //germline
{
	//calculate variant statistics (after annotation because it needs the ID and ANN fields)
	$parser->exec(get_path("ngs-bits")."VariantQC", "-in $annfile -out $stafile", true);
	
	$args = array("-in $annfile", "-out $varfile");
	if ($multi)
	{
		$args[] = "-multi";
	}
	$parser->execTool("NGS/vcf2gsvar.php", implode(" ", $args));
}
else //somatic
{
	$extra = "-t_col $t_col";
	if($n_col!="na" && !empty($n_col))	$extra .= " -n_col $n_col";
	
	//annotate additional columns, somatic only
	$tmp_annfile = $parser->tempFile("_somatic.vcf");
	if(!rename($annfile,$tmp_annfile))	trigger_error("Could not move '$annfile' to '$tmp_annfile'.",E_USER_ERROR);
	$cols = array("Interpro_domain");
	$parser->exec(get_path("SnpSift"), "dbnsfp -noLog -db ".get_path("data_folder")."/dbs/dbNSFP/dbNSFPv2.9.3.txt.gz -f ".implode(",",$cols)." $tmp_annfile > $annfile", true);
	//SnpSift vcf comments are printed twice using dbNSFP -> remove duplicate comments
	$af = Matrix::fromTSV($annfile);
	$af->setComments(array_unique($af->getComments()));
	$af->toTSV($annfile);
}

//zip annotated VCF file
$parser->exec("bgzip", "-c $annfile > $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-p vcf $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

//use exonic/splicing variant list for WGS only (otherwise the NGSD annotation takes too long)
if ($sys['type']=="WGS" && ($sys['build']=="hg19" || $sys['build']=="GRCh37")) 
{
	rename($varfile, $varfile_full);
	$tmp = $parser->tempFile(".bed");
	file_put_contents($tmp, "chrMT\t0\t16569");
	$roi_with_mito = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedAdd", "-in ".get_path("data_folder")."/enrichment/ssHAEv6_2017_01_05.bed -in2 {$tmp} -out {$roi_with_mito}", false);
	$parser->exec(get_path("ngs-bits")."BedMerge", "-in {$roi_with_mito} -out {$roi_with_mito}", false);
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $varfile_full -out $varfile -reg {$roi_with_mito}", true);
}

//annotated variant frequencies from NGSD (not for somatic)
if(!$no_ngsd && $t_col=="na")
{
	//find processed sample with equal processing system for NGSD-annotation
	$extras = array();
	$extras[] = "-psname $out_name";
	if(get_processed_sample_id($out_name, false)==-1)
	{
		$extras = array();
		$tmp = get_processed_sample_name_by_processing_system($sys['name_short'], false);
		if($tmp !== false)
		{
			$extras[] = "-psname $tmp";
			trigger_error("No valid processed sample id found ($out_name)! Used processing id $tmp instead. Statistics may be skewed!", E_USER_WARNING);
		}
	}
	$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $varfile -out $varfile ".implode(" ", $extras), true);
}

//check output TSV file (not for somatic)
if(!$no_fc && $t_col=="na")
{
	$parser->execTool("NGS/check_tsv.php", "-in $varfile");
}

?>
