<?php

/**
  @page vc_sawfish
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_sawfish", "Call of structural variants with sawfish. Creates an VCF file.");
$parser->addInfileArray("bams", "Indexed and sorted BAM file(s).", false);
$parser->addStringArray("genders", "Genders of the provided samples.", false);
$parser->addOutfile("out", "Output VCF file (gzipped and tabix indexed).", false);

//optional
$parser->addStringArray("sample_ids", "Optional sample id(s)/name(s) for VCF output.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addInfileArray("small_var_vcfs", "Optional list of small variant VCFs for each sample.", true);
$parser->addFlag("keep_common_cnvs", "Do not filter out common CNVs (only GRCh38).");
$parser->addInt("threads", "Number of threads used.", true, 4);
$parser->addFlag("test", "Run tool in test-mode (e.g. disable CNV calling)");

extract($parser->parse($argv));

//get sample names
if(!is_null($sample_ids) && (count($sample_ids) != 0))
{
	//check provided sample names
	if(count($sample_ids) != count($bams)) trigger_error("Number of provided BAM files and sample names does not match!", E_USER_ERROR);
}
else
{
	$sample_ids = array();
	foreach($bams as $single_sample_bam)
	{
		$sample_ids[] = basename2($single_sample_bam);
	}
}

//check genders
if(count($bams) != count($genders)) trigger_error("Number of provided BAM files and genders does not match!", E_USER_ERROR);
foreach ($genders as $gender) 
{
	$valid_genders = array("male", "female", "n/a");
	if (array_search($gender, $valid_genders) === false) trigger_error("List of genders contain invalid gender '$gender'!", E_USER_ERROR);
}

//check small variant VCFs
if(!is_null($small_var_vcfs) && (count($small_var_vcfs) != 0))
{
	//check provided VCFs
	if(count($small_var_vcfs) != count($bams)) trigger_error("Number of provided BAM and small variant VCF files does not match!", E_USER_ERROR);
}

$genome = genome_fasta($build);


$discover_output = [];
$jointcall_output = $parser->tempFolder("_sawfish_joint_call");;

//run discover step
for ($i=0; $i < count($bams); $i++) 
{ 
	$discover_output[] = $parser->tempFolder("_sawfish_discover");
	$args = [];
	$in_files = [];
	$args[] = "--threads ".$threads;
	$args[] = "--ref ".$genome;
	$in_files[] = $build;
	$args[] = "--bam ".$bams[$i];
	$in_files[] = $bams[$i];
	if ($build == "GRCh38") 
	{
		if ($genders[$i] == "female") $args[] = "--expected-cn /opt/sawfish-v2.2.1-x86_64-unknown-linux-gnu/data/expected_cn/expected_cn.hg38.XX.bed";
    	if ($genders[$i] == "male") $args[] = "--expected-cn /opt/sawfish-v2.2.1-x86_64-unknown-linux-gnu/data/expected_cn/expected_cn.hg38.XY.bed";
		if ($keep_common_cnvs) $args[] = "--cnv-excluded-regions /opt/sawfish-v2.2.1-x86_64-unknown-linux-gnu/data/cnv_excluded_regions/annotation_only.hg38.bed.gz";
		else $args[] = "--cnv-excluded-regions /opt/sawfish-v2.2.1-x86_64-unknown-linux-gnu/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz";
	}
	if ($test) $args[] = "--disable-cnv"; //disable CNV calling for testing (will not work on small testset)
	if(!is_null($small_var_vcfs) && (count($small_var_vcfs) != 0)) $args[] = "--maf ".$small_var_vcfs[$i];
	$args[] = "--output-dir ".$discover_output[$i];
	$parser->execApptainer("sawfish", "sawfish discover", implode(" ", $args), $in_files);
}

//run joint-call step
$args = [];
$in_files = [];
$args[] = "--threads ".$threads;
$args[] = "--ref ".$genome;
$in_files[] = $build;
for ($i=0; $i < count($bams); $i++) 
{
	$args[] = "--sample ".$discover_output[$i];
	$in_files[] = $discover_output[$i];
	$in_files[] = $bams[$i];
	if(!is_null($small_var_vcfs) && (count($small_var_vcfs) != 0)) $in_files[] = $small_var_vcfs[$i];
}
$args[] = "--output-dir ".$jointcall_output;
$parser->execApptainer("sawfish", "sawfish joint-call", implode(" ", $args), $in_files);

//TODO: Post-processing?

//copy final VCF to output
$parser->copyFile($jointcall_output."/genotyped.sv.vcf.gz", $out);
$parser->copyFile($jointcall_output."/genotyped.sv.vcf.gz.tbi", $out.".tbi");

?>