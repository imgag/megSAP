<?php 
/** 
	@page vc_clair
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_clair_mosaic", "Mosaic variant calling with Clair-Mosaic.");
$parser->addInfile("bam",  "Input file in BAM/CRAM format. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output VCF.GZ file.", false);
$parser->addString("name", "Base file name, typically the processed name ID (e.g. 'GS120001_01').", false);
$parser->addString("model", "clair-mosaic model name used for calling.", false);

//optional
$parser->addInfile("target",  "Enrichment targets BED file.", true);
$parser->addInt("target_extend",  "Call variants up to n bases outside the target region (they are flagged as 'off-target' in the filter column).", true, 0);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addFlag("keep_temp_data", "Keep clair temp folder in sample folder.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");

extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//set bind paths for clair-mosaic container
$in_files = array();
$out_files = array();


//output files
if ($keep_temp_data)
{
	$clair_temp = dirname($out)."/clair_mosaic_temp";
	$out_files[] = $clair_temp;
}
else
{
	$clair_temp = $parser->tempFolder("clair_mosaic_temp");
}


//create basic variant calls
$args = array();
$args[] = "--bam_fn {$bam}";
$args[] = "--ref_fn $genome";
$args[] = "--output_dir {$clair_temp}";
$args[] = "--threads={$threads}";
$args[] = "--platform {$model}";

$args[] = "--sample_name {$name}";
$args[] = "--enable_indel_calling";

// add input files
$in_files[] = $bam;
$in_files[] = $genome;



//calculate target region
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
		$parser->execApptainer("ngs-bits", "BedAdd", "-in $target_extended ".repository_basedir()."/data/misc/special_regions.bed -out $target_extended", [repository_basedir()."/data/misc/special_regions.bed"]);
	}
	
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->execApptainer("ngs-bits", "BedMerge", "-in $target_extended -out $target_merged");
	
	$args[] = "--bed_fn {$target_merged}";
}

//run Clair3 container
$parser->execApptainer("clair-mosaic", "/opt/bin/run_clair_mosaic", implode(" ", $args), $in_files, $out_files, false, true, true, true);
$clair_vcf = $clair_temp."/merge_output.vcf.gz";

//merge SNVs/Indels
$snvs = $clair_temp."/snv.vcf.gz";
$indels = $clair_temp."/indel.vcf.gz";
$clair_merged_vcf1 = $parser->tempFile("_merged.vcf");
$clair_merged_vcf2 = $parser->tempFile("_merged_sorted.vcf.gz");
$parser->execApptainer("ngs-bits", "VcfAdd", "-in {$snvs} {$indels} -out {$clair_merged_vcf1}", [$snvs, $indels]);
$parser->execApptainer("ngs-bits", "VcfSort", "-compression_level 5 -in {$clair_merged_vcf1} -out {$clair_merged_vcf2}");

//post-processing 
$pipeline = array();

//stream vcf
$pipeline[] = array("zcat", "{$clair_merged_vcf2}");

//filter variants according to variant quality>5
$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfFilter", "-remove_invalid -ref $genome", [$genome], [], true));

//split complex variants to primitives
//this step has to be performed before VcfBreakMulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = ["", $parser->execApptainer("vcflib", "vcfallelicprimitives", "-kg", [], [], true)];

//split multi allelic variants
$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "", [], [], true));
//normalize all variants and align INDELs to the left
$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref $genome", [$genome], [], true));

//sort variants by genomic position
$uncompressed_vcf = $parser->tempFile(".vcf");
$pipeline[] = array("", $parser->execApptainer("ngs-bits", "VcfStreamSort", "", [], [], true));

//fix error in VCF file and strip unneeded information
$pipeline[] = array("php ".repository_basedir()."/src/Tools/vcf_fix.php", "--clair3_mode > {$uncompressed_vcf}", false);

//execute post-processing pipeline
$parser->execPipeline($pipeline, "clair post processing");

//add name/pipeline info to VCF header
$vcf = Matrix::fromTSV($uncompressed_vcf);
$comments = $vcf->getComments();
$comments[] = "#reference={$genome}\n";
$comments[] = "#fileDate=".date("Ymd")."\n";
$comments[] = "#ANALYSISTYPE=GERMLINE_SINGLESAMPLE\n";
$comments[] = "#PIPELINE=".repository_revision(true)."\n";
$comments[] = gsvar_sample_header($name, array("DiseaseStatus"=>"Affected"), "#", "");
$vcf->setComments($comments);
$vcf->toTSV($uncompressed_vcf);

//mark off-target variants
if ($target_extend>0)
{
	$on_target_region = $parser->tempFile(".bed");
	$result = $parser->execApptainer("ngs-bits", "BedAdd", "-in {$target} {$target_mito} -out {$on_target_region}", [$target]);
	$tmp = $parser->tempFile(".vcf");
	$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in $uncompressed_vcf -mark off-target -reg {$on_target_region} -out $tmp", [$on_target_region]);
	$parser->execApptainer("htslib", "bgzip", "-c $tmp > $out", [], [dirname($out)]);
}
else
{
	$parser->execApptainer("htslib", "bgzip", "-c $uncompressed_vcf > $out", [], [dirname($out)]);
}

//index output file
$parser->execApptainer("htslib", "tabix", "-f -p vcf $out", [], [dirname($out)]);

?>
