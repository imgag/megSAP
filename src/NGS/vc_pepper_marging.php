<?php 
/** 
	@page vc_pepper_marging
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_pepper_marging", "Variant calling with PepperMargin and DeepVariant.");
$parser->addInfile("bam",  "Input files in BAM format. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output file in VCF.GZ format.", false);
//optional
$parser->addFlag("enable_phasing", "Output phased variants using whatshap.");
$parser->addInfile("target",  "Enrichment targets BED file.", true);
$parser->addInt("target_extend",  "Call variants up to n bases outside the target region (they are flagged as 'off-target' in the filter column).", true, 0);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
// $parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.15);
// $parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 1);
// $parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 15);
// $parser->addInt("min_ao", "Minimum alternative base observation count.", true, 3);
// $parser->addFlag("no_ploidy", "Use freebayes parameter -K, i.e. output all alleles which pass input filters, regardles of genotyping outcome or model.");
// $parser->addFlag("no_bias", "Use freebayes parameter -V, i.e. ignore strand bias and read end distance bias.");
// $parser->addInt("min_qsum", "Minimum quality sum used for variant calling.", true, 0);
// $parser->addFlag("raw_output", "return the raw output of freebayes with no post-processing.");
extract($parser->parse($argv));


//TODO: run as Docker

//init
$genome = genome_fasta($build);
$out_folder = $parser->tempFolder("clair3");

//TODO: implement model path
$model_path = "TODO";

//create basic variant calls
$args = array();
$args[] = "--bam_fn={$bam}";
$args[] = "--ref_fn={$genome}";
$args[] = "--model_path={$model_path}";
$args[] = "--threads={$threads}";
$args[] = "--platform=\"ont\"";
$args[] = "--output={$out_folder}";
if($enable_phasing) $args[] = "--enable_phasing";

//calculate target region
if(isset($target))
{	
	//extend by 'n' bases
	$target_extended = $parser->tempFile("_roi_extended.bed");
	if ($target_extend>0)
	{
		$parser->exec(get_path("ngs-bits")."BedExtend"," -in $target -n $target_extend -out $target_extended -fai {$genome}.fai", true);
	}
	else
	{
		$parser->copyFile($target, $target_extended);
	}
	
	//add special target regions (regions with known pathogenic variants that are often captured by exome/panel, but not inside the target region)
	if ($build=="GRCh38" && $target_extend>0) //only if extended (otherwise it is also added for chrMT calling, etc.)
	{
		$parser->exec(get_path("ngs-bits")."BedAdd"," -in $target_extended ".repository_basedir()."/data/misc/special_regions.bed -out $target_extended ", true);
	}
	
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->exec(get_path("ngs-bits")."BedMerge"," -in $target_extended -out $target_merged", true);
	
	$args[] = "--bed_fn={$target_merged}";
}

//run Clair3
$parser->exec(get_path("clair3")."/run_clair3.sh", implode(" ", $args));


//post-processing 
$pipeline = array();

//TODO: stream vcf.gz
$pipeline[] = array("zcat", "{$vcf_file}");

//filter variants according to variant quality>5
$pipeline[] = array(get_path("vcflib")."vcffilter", "-f \"QUAL > 5\"");

//split complex variants to primitives
//this step has to be performed before vcfbreakmulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = array(get_path("vcflib")."vcfallelicprimitives", "-kg");

//split multi-allelic variants
$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");

//normalize all variants and align INDELs to the left
$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize", "-stream -ref $genome");

//sort variants by genomic position
$pipeline[] = array(get_path("ngs-bits")."VcfStreamSort", "");

//fix error in VCF file and strip unneeded information
$pipeline[] = array("php ".repository_basedir()."/src/NGS/vcf_fix.php", "", false);

//zip
$pipeline[] = array("bgzip", "-c > $out", false);

//(2) execute pipeline
$parser->execPipeline($pipeline, "freebayes post processing");

//(3) mark off-target variants
if ($target_extend>0)
{
	$tmp = $parser->tempFile(".vcf");
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $out -mark off-target -reg $target -out $tmp", true);
	$parser->exec("bgzip", "-c $tmp > $out", false);
}

//(4) index output file
$parser->exec("tabix", "-f -p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return

?>
