<?php 
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_varscan2", "Variant calling with VarScan2.");
$parser->addInfile("bam",  "Input file in BAM format. Space separated. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output file in VCF.GZ format.", false);
//optional
$parser->addInfile("target",  "Enrichment target BED file.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.05);
$parser->addInt("min_dp", "Minimum depth cutoff for variant calling.", true, 20);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 15);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 1);
$parser->addFloat("pval_thres", "p-Value threshold for varinat calling", true, 0.0001);
$parser->addString("name", "Sample name to be used in output.", true, "SAMPLE");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//Call SNPs
$tmp_snp_file = temp_file("varscan2_snps.vcf");
$pipeline = array(
	array(get_path("samtools"), "mpileup --min-MQ $min_mq -f $genome $bam"),
	array(get_path("varscan2"), "mpileup2snp --min-var-freq $min_af --min_avg_qual $min_bq --min-reads2 3 --min_coverage $min_dp --output-vcf --p-value $pval_thres"),
	array(get_path("ngs-bits")."/VcfLeftNormalize", "-ref $genome -out $tmp_snp_file")
);

$parser->execPipeline($pipeline, "Varscan2 SNPs");

//Call INDELs
$tmp_indel_file = temp_file("varscan2_indels.vcf");
$pipeline = array(
	array(get_path("samtools"), "mpileup --min-MQ $min_mq -f $genome $bam" ),
	array(get_path("varscan2"), "mpileup2indel --min-var-freq $min_af --min_avg_qual $min_bq --min-reads2 3 --min_coverage $min_dp --output-vcf --p-value $pval_thres"),
	array(get_path("ngs-bits")."/VcfLeftNormalize", "-ref $genome"),
	array("egrep", "-v '##|#CHROM' > $tmp_indel_file") //filter out header lines (already included in calls for snps, will be merged in later step)
);
$parser->execPipeline($pipeline, "Varscan2 INDELs");

//merge/rename header SNPs and INDELs into one vcf file
$tmp_merged_file = temp_file("varscan2_all_variants.vcf");
$pipeline = array(
	array("cat", "$tmp_snp_file $tmp_indel_file"),
	array("sed", "'s/Sample1/{$name}/' > $tmp_merged_file"), //replace sample1 by actual sample $name
);
$parser->execPipeline($pipeline, "merge");

//Sort file
$tmp_sorted_file = temp_file("varscan2_all_variants_sorted.vcf");
$parser->exec(get_path("ngs-bits")."/VcfSort", "-in $tmp_merged_file -out $tmp_sorted_file");
$parser->exec("bgzip", "-c $tmp_sorted_file > $out");


//flag off-targets
if(isset($target))
{
	$tmp = temp_file("target.vcf");
	$parser->exec(get_path("ngs-bits") . "/VariantFilterRegions", "-in $out -reg $target -mark off-target -out $tmp");
	$parser->exec("bgzip", "-c $tmp > $out");
}

//index
$parser->exec("tabix", "-p vcf $out", false); 

?>
