<?php
/**
  @page vc_delly
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_delly", "Call germline structural variants with delly. Creates a VCF file.");
$parser->addOutfile("out","Output vcf.gz file",false);
$parser->addInFile("t_bam","Tumor BAM file for somatic mode.",true);
$parser->addInfileArray("bam", "Normal BAM file(s). Only one normal BAM file allowed for somatic mode.", false, false);
$parser->addInt("max_threads","Maximum number of threads to use. Should not exceed the number of samples.",true,2);
$parser->addString("build","The genome build to use.",true,"GRCh37");
$parser->addInFile("target","Enrichment target .bed file.",true);
$parser->addInFile("exclude","Exclude of e.g. telomeric or centromeric regions",true);
extract($parser->parse($argv));

//Max number of threads according max number of samples
$threads = (count($bam) < $max_threads) ? count($bam) : $max_threads;  


//Reference genome
$genome = genome_fasta($build);

//switch for somatic and germline mode
$somatic = isset($t_bam);

//original bcf out file of delly
$tmp_bcf_out = $parser->tempFile(".bcf");

if(!$somatic) //germline SV calling
{
	//Set openMP env variable for max threads
	putenv("OMP_NUM_THREADS={$threads}");
	
	//Execute Delly
	$args = [
	"-g",$genome,
	"-o",$tmp_bcf_out
	];
	if(isset($exclude)) $args[] = "-x $exclude";
	$args[] = 	implode(" ",$bam);
	
	$parser->exec(get_path("delly"). " call",implode(" ",$args),true,true);
}
else //somatic SV calling
{
	if(count($bam) != 1)
	{
		trigger_error("More than one normal file specified in tumor-normal mode. Aborting.",E_USER_ERROR);
	}
	putenv("OMP_NUM_THREADS={$threads}");
	
	//Execute Delly
	$args = [
	"-g",$genome,
	"-o",$tmp_bcf_out,
	];
	
	if(isset($exclude)) $args[] = "-x $exclude";
	
	$args[] = $t_bam;
	$args[] = implode(" ",$bam);
	
	$parser->exec(get_path("delly"). " call",implode(" ",$args),true,true);
}

/********************************
 * PARSE RAW BCF FILE TO VCF.GZ *
 ********************************/
//transform to VCF
$tmp_vcf_file = $parser->tempFile("SV_raw.vcf");
$parser->exec(get_path("bcftools")." view"," $tmp_bcf_out > $tmp_vcf_file",true,true);

//sort VCF
$tmp_vcf_sorted = $parser->tempFile("SV_sorted.vcf");
$parser->exec(get_path("ngs-bits")."VcfSort","-in $tmp_vcf_file -out $tmp_vcf_sorted", true,true);

//Flag outliers out of target region
$tmp_vcf_filtered = $parser->tempFile("SV_filtered.vcf");
if(isset($target))
{
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $tmp_vcf_sorted -mark off-target -reg $target -out $tmp_vcf_filtered", true);
}
else
{
	$tmp_vcf_filtered = $tmp_vcf_sorted;
}

//Compress and index
$parser->exec("bgzip -c $tmp_vcf_filtered > $out",true,true);
$parser->exec("tabix", "-p vcf $out", true,true);
?>
