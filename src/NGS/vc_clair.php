<?php 
/** 
	@page vc_clair
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_clair", "Variant calling with Clair3.");
$parser->addInfile("bam",  "Input files in BAM format. Note: .bam.bai file is required!", false);
$parser->addString("folder", "Destination folder for output files.", false);
$parser->addString("name", "Base file name, typically the processed name ID (e.g. 'GS120001_01').", false);
$parser->addInfile("target",  "Enrichment targets BED file.", false);
$parser->addInfile("model", "Model file used for calling.", false);

//optional
$parser->addInt("target_extend",  "Call variants up to n bases outside the target region (they are flagged as 'off-target' in the filter column).", true, 0);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");

extract($parser->parse($argv));

//TODO: check output folder

//init
$genome = genome_fasta($build);

//output files
$clair_temp = "{$folder}/clair_temp";
$out = "{$folder}/{$name}_var.vcf.gz";
$out_gvcf = "{$folder}/{$name}_var.gvcf.gz";

//create basic variant calls
$args = array();
$args[] = "--bam_fn={$bam}";
$args[] = "--ref_fn={$genome}";
$args[] = "--model_path={$model}";
$args[] = "--threads={$threads}";
$args[] = "--platform=\"ont\"";
$args[] = "--keep_iupac_bases";
$args[] = "--gvcf";
//TODO: move to temp
$args[] = "--output={$clair_temp}";

//define env parameter
$args[] = "--sample_name={$name}";
$args[] = "--samtools=".get_path("samtools");
$args[] = "--python=".get_path("python3");
$args[] = "--pypy=".get_path("pypy3");
$args[] = "--parallel=".get_path("parallel");
$args[] = "--whatshap=".get_path("whatshap");

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
putenv("PYTHONPATH=".dirname(get_path("clair3")));
$parser->exec(get_path("clair3"), implode(" ", $args));
$clair_vcf = $clair_temp."/merge_output.vcf.gz";
$clair_gvcf = $clair_temp."/merge_output.gvcf.gz";

//post-processing 
$pipeline = array();

//stream vcf.gz
$pipeline[] = array("zcat", "{$clair_vcf}");

//filter variants according to variant quality>5
$pipeline[] = array(get_path("ngs-bits")."VcfFilter", "-qual 5 -remove_invalid");

//split complex variants to primitives
//this step has to be performed before vcfbreakmulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = array(get_path("vcflib")."vcfallelicprimitives", "-kg");

//split multi-allelic variants
$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");

//normalize all variants and align INDELs to the left
$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize", "-stream -ref $genome");

//sort variants by genomic position
$uncompressed_vcf = $parser->tempFile(".vcf");
$pipeline[] = array(get_path("ngs-bits")."VcfStreamSort", "");

//fix error in VCF file and strip unneeded information
$pipeline[] = array("php ".repository_basedir()."/src/NGS/vcf_fix.php", "--longread_mode > {$uncompressed_vcf}", false);

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
$vcf->setComments(sort_vcf_comments($comments));
$vcf->toTSV($uncompressed_vcf);

//mark off-target variants
if ($target_extend>0)
{
	$tmp = $parser->tempFile(".vcf");
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $uncompressed_vcf -mark off-target -reg $target -out $tmp", true);
	$parser->exec("bgzip", "-c $tmp > $out", false);
}
else
{
	$parser->exec("bgzip", "-c $uncompressed_vcf > $out", false);
}

//index output file
$parser->exec("tabix", "-f -p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return


//create/copy gvcf:
$pipeline = array();
$pipeline[] = array("zcat", $clair_gvcf);
$pipeline[] = array(get_path("ngs-bits")."VcfFilter", "-remove_invalid");
$pipeline[] = array("bgzip", "-c > {$out_gvcf}");
$parser->execPipeline($pipeline, "gVCF post processing");

$parser->exec("tabix", "-f -p vcf {$out_gvcf}", false); //no output logging, because Toolbase::extractVersion() does not return

?>
