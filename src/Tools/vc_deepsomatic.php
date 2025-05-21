<?php 
/** 
	@page vc_deepsomatic
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_deepsomatic", "Variant calling with DeepVariant.");
$parser->addInfile("bam_tumor",  "Input tumor reads. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output file in VCF.gz format.", false);
$parser->addString("model_type", "Type of model to use for variant calling. Choose from <WGS|WES|PACBIO|ONT|FFPE_WGS|FFPE_WES|WGS_TUMOR_ONLY|PACBIO_TUMOR_ONLY|ONT_TUMOR_ONLY>.", false);
//optional
$parser->addInfile("bam_normal", "Input normal reads. Note: .bam.bai file is required!", true);
$parser->addInfile("target",  "Enrichment targets BED file.", true);
$parser->addString("normal_id", "Sample name for normal sample. If not specified, will be inferred from the header information from --bam_tumor", true);
$parser->addString("tumor_id", "Sample name for tumor sample. If not specified, will be inferred from the header information from --bam_tumor", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("gvcf", "Enable output of gVCF files and define output filepath for gVCF files.", true, "");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.15);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 20);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 10);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addFlag("raw_output", "return the raw output of deepvariant with no post-processing.");
$parser->addFlag("allow_empty_examples", "allows DeepVariant to call variants even if no examples were created with make_examples.");
$parser->addFlag("tumor_only", "run DeepSomatic in tumor-only mode");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);
if (empty($tumor_id)) $tumor_id = basename2($bam_tumor);
if (!$tumor_only && empty($normal_id)) $normal_id = basename2($bam_normal);

//create basic variant calls
$args = array();
$in_files = array();

$args[] = "--regions $target";
$in_files[] = $target; 

if ($allow_empty_examples)
{
	$args[] = "--call_variants_extra_args=allow_empty_examples=true";
}

$args[] = "--model_type=$model_type";
$args[] = "--make_examples_extra_args=min_mapping_quality=$min_mq,min_base_quality=$min_bq,vsc_min_fraction_indels=$min_af,vsc_min_fraction_snps=$min_af";
$args[] = "--ref=$genome";
$args[] = "--reads_tumor=$bam_tumor";
$args[] = "--sample_name_tumor=$tumor_id";
$args[] = "--num_shards=".$threads;

if (!empty($gvcf)) $args[] = "--output_gvcf=$gvcf";
if (!$tumor_only) 
{
	$args[] = "--reads_normal=$bam_normal";
	$args[] = "--sample_name_normal=$normal_id";
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

if (!$tumor_only)
{
	//################################################################################################
	//Filter variants
	//################################################################################################

	$variants = Matrix::fromTSV($tmp_out);
	$variants_filtered = new Matrix();

	//fix column names
	$colnames = $variants->getHeaders();
	$colidx_tumor = array_search($tumor_id, $colnames);
	$colidx_normal = array_search($normal_id, $colnames);

	//quality cutoffs (taken from Strelka)
	$min_td = 20;
	$min_taf = 0.05;
	$min_tsupp = 3;
	$min_nd = 15;
	$max_naf_rel = 1/6;
	$min_call_sq = 5;
	$min_filter_sq = 17.5;

	//set comments and column names
	$filter_format = '#FILTER=<ID=%s,Description="%s">';
	$comments = [
		sprintf($filter_format, "all-unknown", "Allele unknown"),
		sprintf($filter_format, "special-chromosome", "Special chromosome"),
		sprintf($filter_format, "depth-tum", "Sequencing depth in tumor is too low (< {$min_td})"),
		sprintf($filter_format, "freq-tum", "Allele frequency in tumor < {$min_taf}"),
		sprintf($filter_format, "lt-3-reads", "Less than {$min_tsupp} supporting tumor reads")
		];

	$variants_filtered->setComments(sort_vcf_comments(array_merge($variants->getComments(), $comments)));
	$variants_filtered->setHeaders($colnames);

	for($i = 0; $i < $variants->rows(); ++$i)
	{
		$row = $variants->getRow($i);

		$ref = $row[3];
		$alt = $row[4];
		$format = $row[8];
		$tumor = $row[$colidx_tumor];
		$normal = $row[$colidx_normal];

		$filters = [];
		$type = (strlen($row[3]) > 1 || strlen($row[4]) > 1) ? "INDEL" : "SNV";

		$filter = array_diff(explode(";", $row[6]), ["PASS"]);
		//filter out variants called as germline variants
		if (in_array("GERMLINE", $filter)) continue;

		if (!preg_match("/^[acgtACGT]*$/", $alt))
		{
			$filter[] = "all-unknown";
		}
		if (chr_check($row[0], 22, false) === FALSE)
		{
			$filter[] = "special-chromosome";
		}
		$calls = [];
		list($td, $tf) = vcf_deepvariant($format, $tumor);
		$calls[] = [ $alt, $td, $tf, $filter];

		foreach ($calls as $call)
		{
			$variant = $row;
			list($alt, $td, $tf, $filter) = $call;
			$variant[4] = $alt;

			if ($td * $tf < $min_tsupp) $filter[] = "lt-3-reads";
			if ($td < $min_td) $filter[] = "depth-tum";
			if ($tf < $min_taf) $filter[] = "freq-tum";

			if (empty($filter))
			{
				$filter[] = "PASS";
			}

			$variant[6] = implode(";", $filter);
			$variants_filtered->addRow($variant);		
		}
	}
	$final = $parser->tempFile("_filtered.vcf");
	$variants_filtered->toTSV($final);
}
else $final = $tmp_out;

//flag off-target variants
if (!empty($target))
{
	if (!$tumor_only)
	{
		//add mito region to target region (for WGS):
		$target_mito = $parser->tempFile("_mito.bed");
		file_put_contents($target_mito, "chrMT\t0\t16569");
		$full_target = $parser->tempFile("_full.bed");
		$parser->execApptainer("ngs-bits", "BedAdd", "-in {$target} {$target_mito} -out {$full_target}", [$target]);
		$target = $full_target;
	}

	$vcf_offtarget = $parser->tempFile("_filtered.vcf");
	$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in $final -mark off-target -reg $target -out $vcf_offtarget", [$target]);
	$final = $vcf_offtarget;
}

//remove artefacts specific for processing system (blacklist)
//artefacts are caused e.g. by hairpin sequences when using enzymatic digestion
//how artefacts are determined is documented in /mnt/storage3/users/ahsturm1/Sandbox/2023_01_31_twist_indel_artefacts/
if (!empty($target) && !$tumor_only)
{
	$artefact_vcf = repository_basedir()."/data/misc/enzymatic_digestion_artefacts/".basename($target, ".bed").".vcf";
	if (file_exists($artefact_vcf))
	{
		$vcf_with_artefacts = dirname($out)."/".basename($out, ".vcf.gz")."_with_enzymatic_artefacts.vcf.gz";
		$parser->execApptainer("htslib", "bgzip", "-c $final > $vcf_with_artefacts", [], [dirname($vcf_with_artefacts)]);
		$vcf_no_artefacts = $parser->tempFile("_filtered_no_artefacts.vcf");
		$parser->execApptainer("ngs-bits", "VcfSubtract", "-in $final -in2 $artefact_vcf -out $vcf_no_artefacts", [$artefact_vcf]);
		$final = $vcf_no_artefacts;
	}
}

//zip
$parser->execApptainer("htslib", "bgzip", "-c $final > $out", [], [dirname($out)]);

//(3) index output file
$parser->execApptainer("htslib", "tabix", "-f -p vcf $out", [], [dirname($out)]);

?>
