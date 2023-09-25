<?php

/**
  @page vc_sniffles
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_sniffles", "Call of structural variants with sniffles. Creates an VCF file.");
$parser->addInfile("bam", "Indexed and sorted BAM file.", false);
$parser->addOutfile("out", "Output VCF file (gzipped and tabix indexed).", false);
//optional
$parser->addString("name", "Optional sample name for VCF output.", true, "");
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFlag("somatic", "Use somatic mode for SV calling.");
$parser->addInt("threads", "Number of threads used.", true, 4);
extract($parser->parse($argv));

//TODO: add tandem repeat file

if($name == "") $name = basename2($bam);
$tmp_vcf = $parser->tempFile(".vcf", "sniffles");

$args = array();
$args[] = "--input ".$bam;
$args[] = "--vcf ".$tmp_vcf;
$args[] = "--threads ".$threads;
$args[] = "--reference ".genome_fasta($build);
$args[] = "--allow-overwrite";
// $args[] = "--output-rnames ".$name;
if($somatic) $args[] = "--non-germline";

$sniffles = get_path("sniffles");

$parser->exec(get_path("sniffles"), implode(" ", $args), true);

//add name/pipeline info to VCF header
$vcf = Matrix::fromTSV($tmp_vcf);
$comments = $vcf->getComments();
$comments[] = "#reference=".genome_fasta($build)."\n";
$comments[] = "#PIPELINE=".repository_revision(true)."\n";
$comments[] = gsvar_sample_header($name, array("DiseaseStatus"=>"affected"), "#", "");
$vcf->setComments(sort_vcf_comments($comments));
$vcf->toTSV($tmp_vcf);

$parser->exec("bgzip", "-c $tmp_vcf > $out", false);

//index output file
$parser->exec("tabix", "-f -p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return

?>