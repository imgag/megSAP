<?php

/**
	@page annotate
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("annotate", "Annotate variants called on genome build GRCh37.");
$parser->addString("out_name", "Processed sample ID (e.g. 'GS120001_01').", false);
$parser->addString("out_folder", "Output folder.", false);
//optional
$parser->addInfile("system", "Processing system INI file (automatically determined from NGSD if 'out_name' is a valid processed sample name).", true);
$parser->addString("vcf", "Path to (bgzipped) VCF file (if different from {output_folder}/{out_name}_var.vcf.gz).", true, "");
$parser->addFlag("no_fc", "No format check (vcf/tsv).");
$parser->addFlag("multi", "Enable multi-sample mode.");
$parser->addFlag("somatic", "Enable somatic mode (no variant QC and no GSvar file).", true, "na");
$parser->addFlag("updown", "Don't discard up- or downstream anntations (5000 bases around genes).");
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
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
$stafile = $out_folder."/".$out_name."_stats_vc.qcML";

//get system
$sys = load_system($system, $out_name);
if ($sys['build']!="GRCh37")
{
	trigger_error("Unknown genome build ".$sys['build']." cannot be annotated!", E_USER_ERROR);
}

//annotate VCF
$args = array("-in {$vcf_unzipped}", "-out {$annfile}");
$args[] = "-build ".$sys['build'];
$args[] = "-threads {$threads}";
$parser->execTool("NGS/an_vep.php", implode(" ", $args));

//check vcf file
if(!$no_fc)
{
	$parser->exec(get_path("ngs-bits")."VcfCheck", "-in $annfile -ref ".genome_fasta($sys['build']), true);
}

//convert to GSvar file
if (!$somatic) //germline only
{
	//calculate variant statistics (after annotation because it needs the ID and ANN fields)
	$parser->exec(get_path("ngs-bits")."VariantQC", "-in $annfile -out $stafile", true);
	
	$args = array("-in $annfile", "-out $varfile", "-blacklist");
	if ($multi) $args[] = "-genotype_mode multi";
	if ($updown) $args[] = "-updown";
	if ($sys['type']=="WGS") $args[] = "-wgs";
	$parser->execTool("NGS/vcf2gsvar.php", implode(" ", $args));
}

//zip annotated VCF file
$parser->exec("bgzip", "-c $annfile > $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-p vcf $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

//annotated variant frequencies from NGSD (not for somatic)
if(!$somatic && db_is_enabled("NGSD"))
{
	$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in {$varfile} -out {$varfile} -psname {$out_name}", true);
}

//check output TSV file (not for somatic)
if(!$somatic && !$no_fc)
{
	$parser->execTool("NGS/check_tsv.php", "-in $varfile -build ".$sys['build']);
}

?>
