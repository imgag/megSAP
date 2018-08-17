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
$parser->addString("vcf", "Path to (bgzipped) VCF file (if different from {output_folder}/{out_name}_var.vcf.gz).", true, "");
$parser->addFlag("no_fc", "No format check (vcf/tsv).");
$parser->addFlag("multi", "Enable multi-sample mode.");
$parser->addFlag("somatic", "Enable somatic mode (no variant QC and no GSvar file).", true, "na");
$parser->addFlag("updown", "Don't discard up- or downstream anntations (5000 bases around genes).");
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
$varfile_rare = $out_folder."/".$out_name."_rare.GSvar";
$stafile = $out_folder."/".$out_name."_stats_vc.qcML";

//get system
$sys = load_system($system, $out_name);
if ($sys['build']!="hg19" && $sys['build']!="GRCh37" && $sys['build']!="mm10")
{
	trigger_error("Unknown genome build ".$sys['build']." cannot be annotated!", E_USER_ERROR);
}

//annotate VCF
$args = array("-in $vcf_unzipped", "-out $annfile");
$args[] = "-build ".$sys['build'];
$parser->execTool("NGS/an_vep.php", implode(" ", $args));

//check vcf file
if(!$no_fc)
{
	$parser->execTool("NGS/check_vcf.php", "-in $annfile");
}

//convert to GSvar file
if (!$somatic) //germline
{
	//calculate variant statistics (after annotation because it needs the ID and ANN fields)
	$parser->exec(get_path("ngs-bits")."VariantQC", "-in $annfile -out $stafile", true);
	
	$args = array("-in $annfile", "-out $varfile", "-build ".$sys['build']);
	if ($multi) $args[] = "-multi";
	if ($updown) $args[] = "-updown";
	$parser->execTool("NGS/vcf2gsvar.php", implode(" ", $args));
}

//zip annotated VCF file
$parser->exec("bgzip", "-c $annfile > $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-p vcf $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

//use exonic/splicing variant list for WGS only (otherwise the NGSD annotation takes too long)
if ($sys['type']=="WGS" && ($sys['build']=="hg19" || $sys['build']=="GRCh37")) 
{
	$parser->moveFile($varfile, $varfile_full);
	$tmp = $parser->tempFile(".bed");
	file_put_contents($tmp, "chrMT\t0\t16569");
	$roi_with_mito = $parser->tempFile(".bed");
	$parser->exec(get_path("ngs-bits")."BedAdd", "-in ".get_path("data_folder")."/enrichment/ssHAEv6_2017_01_05.bed {$tmp} -out {$roi_with_mito}", false); //TODO use CCDS/Ensembl coding instead! (everywhere)
	$parser->exec(get_path("ngs-bits")."BedMerge", "-in {$roi_with_mito} -out {$roi_with_mito}", false);
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in {$varfile_full} -out {$varfile} -reg {$roi_with_mito}", true);
	$tmp2 = $parser->tempFile(".txt");
	file_put_contents($tmp2, "Allele frequency\tmax_af=1.0");
	$parser->exec(get_path("ngs-bits")."VariantFilterAnnotations", "-in {$varfile_full} -out {$varfile_rare} -filters {$tmp2}", true);
	$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in {$varfile_rare} -out {$varfile_rare} -psname {$out_name}", true);
}

//annotated variant frequencies from NGSD (not for somatic)
if(!$somatic)
{
	$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in {$varfile} -out {$varfile} -psname {$out_name}", true);
}

//check output TSV file (not for somatic)
if(!$somatic && !$no_fc)
{
	$parser->execTool("NGS/check_tsv.php", "-in $varfile");
}

?>
