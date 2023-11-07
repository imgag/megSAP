<?php 
/** 
	@page vc_clair_trio
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_clair_trio", "Trio variant calling with Clair3-Trio.");
$parser->addInfile("bam_c", "BAM file of child (index).", false, true);
$parser->addInfile("bam_f", "BAM file of father.", false, true);
$parser->addInfile("bam_m", "BAM file of mother.", false, true);

$parser->addString("folder", "Destination folder for output files.", false);
$parser->addInfile("target",  "Enrichment targets BED file.", false);
$parser->addInfile("model", "Model file used for calling.", false);
$parser->addInfile("trio_model", "Trio model file used for calling.", false);

//optional
$parser->addInt("target_extend",  "Call variants up to n bases outside the target region (they are flagged as 'off-target' in the filter column).", true, 0);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//output files
$clair_temp = "{$folder}/clair_temp";
$out = "{$folder}/trio_var.vcf.gz";

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
	
	// //add special target regions (regions with known pathogenic variants that are often captured by exome/panel, but not inside the target region)
	// if ($build=="GRCh38" && $target_extend>0) //only if extended (otherwise it is also added for chrMT calling, etc.)
	// {
	// 	$parser->exec(get_path("ngs-bits")."BedAdd"," -in $target_extended ".repository_basedir()."/data/misc/special_regions.bed -out $target_extended ", true);
	// }
	
	// $target_merged = $parser->tempFile("_merged.bed");
	// $parser->exec(get_path("ngs-bits")."BedMerge"," -in $target_extended -out $target_merged", true);
	
	// $args[] = "--bed_fn={$target_merged}";
}

//prepare container

//bind all required paths to the container
$bind_paths = array();
//input dirs (read-only)
$bind_paths[] = dirname($bam_c).":".dirname($bam_c).":ro";
$bind_paths[] = dirname($bam_f).":".dirname($bam_f).":ro";
$bind_paths[] = dirname($bam_m).":".dirname($bam_m).":ro";
$bind_paths[] = dirname($target_extended).":".dirname($target_extended).":ro";
$bind_paths[] = dirname($genome).":".dirname($genome).":ro";
$bind_paths[] = $model.":".$model.":ro";
$bind_paths[] = $trio_model.":".$trio_model.":ro";
//output folder (read-write)
$bind_paths[] = $clair_temp.":".$clair_temp.":rw";
$parser->exec("mkdir", "-p {$clair_temp}");

//set parameters for clair
$args = array();
$args[] = "--bam_fn_c={$bam_c}";
$args[] = "--bam_fn_p1={$bam_f}";
$args[] = "--bam_fn_p2={$bam_m}";
$args[] = "--ref_fn={$genome}";
$args[] = "--model_path_clair3={$model}";
$args[] = "--model_path_clair3_trio={$trio_model}";
$args[] = "--threads={$threads}";
$args[] = "--output={$clair_temp}";
//TODO: test
$args[] = "--gvcf";


//run in container
$parser->execSingularity("clair3-trio", get_path("container_clair3-trio"), $bind_paths, "/opt/bin/run_clair3_trio.sh", implode(" ", $args));

//TODO: check naming
$clair_vcf = $clair_temp."/merge_output.vcf.gz";

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
$comments[] = "#ANALYSISTYPE=GERMLINE_TRIO\n";
$comments[] = "#PIPELINE=".repository_revision(true)."\n";
$comments[] = gsvar_sample_header(basename2($bam_c), array("DiseaseStatus"=>"affected"), "#", "");
$comments[] = gsvar_sample_header(basename2($bam_f), array("DiseaseStatus"=>"unaffected"), "#", "");
$comments[] = gsvar_sample_header(basename2($bam_m), array("DiseaseStatus"=>"unaffected"), "#", "");
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

?>
