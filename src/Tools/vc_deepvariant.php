<?php 
/** 
	@page vc_deepvariant
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_deepvariant", "Variant calling with DeepVariant.");
$parser->addInfile("bam",  "Indexed BAM/CRAM file to call variants on.", false);
$parser->addOutfile("out", "Output file in VCF.GZ format.", false);
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
$parser->addFlag("gpu", "Use DeepVariant container with GPU acceleration support.");
$parser->addFlag("add_sample_header", "Add sample header to VCF file.");
$parser->addString("name", "Sample name for the sample header. Has to be set with add_sample_header.", true, "");
$parser->addString("analysistype", "Type of analysis performed for the sample header. Has to be set with add_sample_header.", true, "GERMLINE_SINGLESAMPLE");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//prepare DeepVariant arguments
$args = [];
$args[] = "--reads={$bam}";
if(isset($target))
{	
	//extend by 'n' bases
	$target_extended = $parser->tempFile("_roi_extended.bed");
	if ($target_extend>0)
	{
		$parser->execApptainer("ngs-bits", "BedExtend", "-in {$target} -n {$target_extend} -out {$target_extended} -fai {$genome}.fai", [$target, $genome]);
	}
	else
	{
		$parser->copyFile($target, $target_extended);
	}
	
	//add special target regions (regions with known pathogenic variants that are often captured by exome/panel, but not inside the target region)
	if ($build=="GRCh38" && $target_extend>0) //only if extended (otherwise it is also added for chrMT calling, etc.)
	{
		$parser->execApptainer("ngs-bits", "BedAdd", "-in {$target_extended} ".repository_basedir()."data/misc/special_regions.bed -out {$target_extended}", [repository_basedir()."data/misc/special_regions.bed"]);
	}
	
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->execApptainer("ngs-bits", "BedMerge", "-in {$target_extended} -out {$target_merged}");
	
	$args[] = "--regions {$target_merged}";
}
if ($allow_empty_examples)
{
	$args[] = "--call_variants_extra_args=allow_empty_examples=true";
}
$args[] = "--model_type={$model_type}";
$args[] = "--make_examples_extra_args=min_mapping_quality={$min_mq},min_base_quality={$min_bq},vsc_min_fraction_indels={$min_af},vsc_min_fraction_snps={$min_af}";
$args[] = "--postprocess_variants_extra_args=cpus={$threads}";
$args[] = "--ref={$genome}";
$args[] = "--num_shards={$threads}";
if (!empty($gvcf))
{
	$args[] = "--output_gvcf={$gvcf}";
}
$args[] = "--intermediate_results_dir=".$parser->tempFolder(); //if not set, examples are written to /tmp/, even if tmp folder is overwritten in environment variables, e.g. in SGE

// run deepvariant
$pipeline = array();
$prefix = container_platform()=='apptainer' ? "APPTAINERENV" : "SINGULARITYENV";
if ($raw_output)
{
	$command = $parser->execApptainer("deepvariant", "run_deepvariant" ,implode(" ", $args)." --output_vcf={$out}", [$genome, $bam], [dirname($out)], true, true, true, true, $gpu);
	$parser->exec("{$prefix}_OMP_NUM_THREADS={$threads} {$prefix}_TF_NUM_INTRAOP_THREADS={$threads} {$prefix}_TF_NUM_INTEROP_THREADS={$threads} {$command}", "");
	return;
}

$vcf_deepvar_out = $parser->tempFile(".vcf.gz");
$command = $parser->execApptainer("deepvariant", "run_deepvariant", implode(" ", $args)." --output_vcf={$vcf_deepvar_out}", [$genome, $bam], [dirname($out)], true, true, true, true, $gpu);
$parser->exec("{$prefix}_OMP_NUM_THREADS={$threads} {$prefix}_TF_NUM_INTRAOP_THREADS={$threads} {$prefix}_TF_NUM_INTEROP_THREADS={$threads} {$command}", "");

//filter variants according to variant quality>5
$pipeline[] = ["zcat", "$vcf_deepvar_out"];
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfFilter", "-qual 5 -remove_invalid -ref $genome", [$genome], [], true)];

//split complex variants to primitives
//this step has to be performed before VcfBreakMulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = ["", $parser->execApptainer("vcflib", "vcfallelicprimitives", "-kg", [], [], true)];

//split multi-allelic variants - -no_errors flag can be removed, when vcfallelicprimitives is replaced
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "-no_errors", [], [], true)];
//normalize all variants and align INDELs to the left
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref $genome", [$genome], [], true)];

//sort variants by genomic position
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfStreamSort", "", [], [], true)];

//fix errors and merge variants
$tmp_out = $parser->tempFile(".vcf");
$pipeline[] = ["php ".repository_basedir()."/src/Tools/vcf_fix.php", "--deepvariant_mode > $tmp_out", false];

//(2) execute pipeline
$parser->execPipeline($pipeline, "deepvariant post processing");

//Add header to VCF file
$vcf = Matrix::fromTSV($tmp_out);
$comments = $vcf->getComments();
$comments[] = "#fileformat=VCFv4.2\n";
$comments[] = "#source=DeepVariant ".get_path("container_deepvariant")."\n";
$comments[] = "#reference=".genome_fasta($build, false)."\n";
$comments[] = "#fileDate=".date("Ymd")."\n";
if ($add_sample_header)
{
	$comments[] = "#ANALYSISTYPE=$analysistype\n";
	$comments[] = "#PIPELINE=".repository_revision(true)."\n";
	$comments[] = gsvar_sample_header($name, array("DiseaseStatus"=>"Affected"), "#", "");
}
$vcf->setComments($comments);
$vcf->toTSV($tmp_out);

//zip
$parser->execApptainer("htslib", "bgzip", "-c $tmp_out > $out", [], [dirname($out)]);

//(3) mark off-target variants
if ($target_extend>0)
{
	$tmp = $parser->tempFile(".vcf");
	$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in $out -mark off-target -reg $target -out $tmp", [$out, $target]);
	$parser->execApptainer("htslib", "bgzip", "-c $tmp > $out", [], [dirname($out)]);
}

//(4) index output file
$parser->execApptainer("htslib", "tabix", "-f -p vcf $out", [], [dirname($out)]);

?>
