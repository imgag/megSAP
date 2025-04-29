<?php 
/** 
	@page vc_deepvariant
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_deepvariant", "Variant calling with DeepVariant.");
$parser->addInfileArray("bam",  "Input files in BAM format. Space separated. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output file in VCF.gz format.", false);
$parser->addString("model_type", "Type of model to use for variant calling. Choose from <WGS|WES|PACBIO|ONT_R104|HYBRID_PACBIO_ILLUMINA|MASSEQ>.", false);
//optional
$parser->addInfile("target",  "Enrichment targets BED file.", true);
$parser->addInt("target_extend",  "Call variants up to n bases outside the target region (they are flagged as 'off-target' in the filter column).", true, 0);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("gvcf", "Enable output of gVCF files and define output filepath for gVCF files.", true, "");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.15);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 20);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 10);
$parser->addFlag("raw_output", "return the raw output of deepvariant with no post-processing.");
$parser->addFlag("allow_empty_examples", "allows DeepVariant to call variants even if no examples were created with make_examples.");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//create basic variant calls
$args = array();
$in_files = array();
if(isset($target))
{	
	//extend by 'n' bases
	$target_extended = $parser->tempFile("_roi_extended.bed");
	if ($target_extend>0)
	{
		$parser->execApptainer("ngs-bits", "BedExtend", "-in $target -n $target_extend -out $target_extended -fai {$genome}.fai", [$target, $genome]);
	}
	else
	{
		$parser->copyFile($target, $target_extended);
	}
	
	//add special target regions (regions with known pathogenic variants that are often captured by exome/panel, but not inside the target region)
	if ($build=="GRCh38" && $target_extend>0) //only if extended (otherwise it is also added for chrMT calling, etc.)
	{
		$parser->execApptainer("ngs-bits", "BedAdd", "-in $target_extended ".repository_basedir()."data/misc/special_regions.bed -out $target_extended ", [repository_basedir()."data/misc/special_regions.bed"]);
	}
	
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->execApptainer("ngs-bits", "BedMerge", "-in $target_extended -out $target_merged");
	
	$args[] = "--regions $target_merged";
}
if ($allow_empty_examples)
{
	$args[] = "--call_variants_extra_args=allow_empty_examples=true";
}

$args[] = "--model_type=$model_type";
$args[] = "--make_examples_extra_args=min_mapping_quality=$min_mq,min_base_quality=$min_bq,vsc_min_fraction_indels=$min_af,vsc_min_fraction_snps=$min_af";
$args[] = "--ref=$genome";
$args[] = "--reads=".implode(" ", $bam);
$args[] = "--num_shards=".$threads;

if (!empty($gvcf)) $args[] = "--output_gvcf=$gvcf";

$in_files[] = $genome;
$in_files = array_merge($in_files, $bam);

// run deepvariant
$pipeline = array();

if ($raw_output)
{
	$parser->execApptainer("deepvariant", "run_deepvariant" ,implode(" ", $args)." --output_vcf=$out", $in_files, [dirname($out)]);
	return;
}

$vcf_deepvar_out = $parser->tempFile(".vcf.gz");
$parser->execApptainer("deepvariant", "run_deepvariant", implode(" ", $args)." --output_vcf=$vcf_deepvar_out", $in_files, [dirname($out)]);

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
$source_line = "##source=DeepVariant ".get_path("container_deepvariant")."\n";
$reference_line = "##reference=".genome_fasta($build, false)."\n";
file_put_contents($tmp_out, $file_format . $file_date . $source_line . $reference_line . file_get_contents($tmp_out));

//zip
$parser->exec("bgzip", "-c $tmp_out > $out");

//(3) mark off-target variants
if ($target_extend>0)
{
	$tmp = $parser->tempFile(".vcf");
	$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in $out -mark off-target -reg $target -out $tmp", [$out, $target]);
	$parser->exec("bgzip", "-c $tmp > $out", false);
}

//(4) index output file
$parser->exec("tabix", "-f -p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return

?>
