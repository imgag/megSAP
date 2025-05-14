<?php 
/** 
	@page vc_deepsomatic
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_deepsomatic", "Variant calling with DeepVariant.");
$parser->addInfileArray("bam_tumor",  "Input tumor reads. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output file in VCF.gz format.", false);
$parser->addString("model_type", "Type of model to use for variant calling. Choose from <WGS|WES|PACBIO|ONT|FFPE_WGS|FFPE_WES|WGS_TUMOR_ONLY|PACBIO_TUMOR_ONLY|ONT_TUMOR_ONLY>.", false);
//optional
$parser->addInfile("bam_normal", "Input normal reads. Note: .bam.bai file is required!", true);
$parser->addInfile("target",  "Enrichment targets BED file.", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("gvcf", "Enable output of gVCF files and define output filepath for gVCF files.", true, "");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.15);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 20);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 10);
$parser->addFlag("raw_output", "return the raw output of deepvariant with no post-processing.");
$parser->addFlag("allow_empty_examples", "allows DeepVariant to call variants even if no examples were created with make_examples.");
$parser->addFlag("tumor_only", "run DeepSomatic in tumor-only mode");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//create basic variant calls
$args = array();
$in_files = array();

$args[] = "--regions $target";
//TODO Kilian add sample names
$args[] = "--sample_name_normal=normal";
$args[] = "--sample_name_tumor=tumor";

if ($allow_empty_examples)
{
	$args[] = "--call_variants_extra_args=allow_empty_examples=true";
}

$args[] = "--model_type=$model_type";
$args[] = "--make_examples_extra_args=min_mapping_quality=$min_mq,min_base_quality=$min_bq,vsc_min_fraction_indels=$min_af,vsc_min_fraction_snps=$min_af";
$args[] = "--ref=$genome";
$args[] = "--reads_tumor=$bam_tumor";
$args[] = "--num_shards=".$threads;

if (!empty($gvcf)) $args[] = "--output_gvcf=$gvcf";
if (!$tumor_only) 
{
	$args[] = "--reads_normal=$bam_normal";
	$in_files[] = $bam_normal;
}

$in_files[] = $genome;
$in_files[] = $bam_tumor; 

// run deepsomatic
$pipeline = array();

if ($raw_output)
{
	$parser->execApptainer("deepsomatic", "run_deepsomatic" ,implode(" ", $args)." --output_vcf=$out", $in_files, [dirname($out)]);
	return;
}

$vcf_deepvar_out = $parser->tempFile(".vcf.gz");
$parser->execApptainer("deepsomatic", "run_deepsomatic", implode(" ", $args)." --output_vcf=$vcf_deepvar_out", $in_files, [dirname($out)]);

//filter variants according to variant quality>5
$pipeline[] = ["zcat", "$vcf_deepvar_out"];
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfFilter", "-qual 5 -remove_invalid -ref $genome", [$genome], [], true)];
/* } */

//split complex variants to primitives
//this step has to be performed before VcfBreakMulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = ["", $parser->execApptainer("vcflib", "vcfallelicprimitives", "-kg", [], [], true)];

//split multi-allelic variants - -no_errors flag can be removed, when vcfallelicprimitives is replaced
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "-no_errors", [], [], true)];
//normalize all variants and align INDELs to the left
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref $genome", [$genome], [], true)];

//sort variants by genomic position
$tmp_out = $parser->tempFile(".vcf");
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfStreamSort", "-out $tmp_out", [], [], true)];

//(2) execute pipeline
$parser->execPipeline($pipeline, "deepvariant post processing");

//prepend source, date and reference to outfile:
$file_format = "##fileformat=VCFv4.2\n";
$file_date = "##fileDate=".date("Ymd")."\n";
$source_line = "##source=DeepSomatic ".get_path("container_deepvariant")."\n";
$reference_line = "##reference=".genome_fasta($build, false)."\n";
file_put_contents($tmp_out, $file_format . $file_date . $source_line . $reference_line . file_get_contents($tmp_out));

//zip
$parser->execApptainer("htslib", "bgzip", "-c $tmp_out > $out", [], [dirname($out)]);

//(3) index output file
$parser->execApptainer("htslib", "tabix", "-f -p vcf $out", [], [dirname($out)]);

?>
