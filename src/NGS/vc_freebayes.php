<?php 
/** 
	@page vc_freebayes
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_freebayes", "Variant calling with freebayes.");
$parser->addInfileArray("bam",  "Input files in BAM format. Space separated. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output file in VCF.GZ format.", false);
//optional
$parser->addInfile("target",  "Enrichment targets BED file.", true);
$parser->addInt("target_extend",  "Call variants up to n bases outside the target region (they are flagged as 'off-target' in the filter column).", true, 0);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.15);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 1);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 10);
$parser->addFlag("no_ploidy", "Use freebayes parameter -K, i.e. output all alleles which pass input filters, regardles of genotyping outcome or model.");
$parser->addFlag("no_bias", "Use freebayes parameter -V, i.e. ignore strand bias and read end distance bias.");
extract($parser->parse($argv));

//(1) set up variant calling pipeline
$genome = get_path("local_data")."/{$build}.fa";
$pipeline = array();

//create basic variant calls
$args = array();
if(isset($target))
{
	if ($target_extend>0)
	{
		$target_extended = $parser->tempFile("_extended.bed");
		$parser->exec(get_path("ngs-bits")."BedExtend"," -in $target -n $target_extend -out $target_extended -fai ".get_path("data_folder")."/genomes/".$build.".fa.fai", true);
	}
	else
	{
		$target_extended = $target;
	}
	
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->exec(get_path("ngs-bits")."BedMerge"," -in $target_extended -out $target_merged", true);
	
	$args[] = "-t $target_merged";
}
if ($no_ploidy)
{
	$args[] = "--pooled-continuous";
}
if ($no_bias)
{
	$args[] = "--binomial-obs-priors-off";
}
$args[] = "--min-alternate-fraction $min_af";
$args[] = "--min-mapping-quality $min_mq";
$args[] = "--min-base-quality $min_bq"; //max 10% error propbability
$args[] = "--min-alternate-qsum 90"; //At least 3 good observations
$pipeline[] = array(get_path("freebayes"), "-b ".implode(" ",$bam)." -f $genome ".implode(" ", $args));

//filter variants according to variant quality>5 , alternate observations>=3
$pipeline[] = array(get_path("vcflib")."vcffilter", "-f \"QUAL > 5 & AO > 2\"");

//split complex variants to primitives
//this step has to be performed before vcfbreakmulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = array(get_path("vcflib")."vcfallelicprimitives", "-kg");

//split multi-allelic variants
$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");

//normalize all variants and align INDELs to the left
$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize","-ref $genome");

//sort variants by genomic position
$pipeline[] = array(get_path("ngs-bits")."VcfStreamSort","");

//fix error in VCF file and strip unneeded information
$pipeline[] = array("php ".repository_basedir()."/src/NGS/vcf_fix.php", "", false);

//zip
$pipeline[] = array("bgzip", "-c > $out", false);

//(2) execute pipeline
$parser->execPipeline($pipeline, "variant calling");

//(3) mark off-target variants
if ($target_extend>0)
{
	$tmp = $parser->tempFile(".vcf");
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $out -mark off-target -reg $target -out $tmp", true);
	$parser->exec("bgzip", "-c $tmp > $out", false);
}

//(4) index output file
$parser->exec("tabix", "-p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return

?>
