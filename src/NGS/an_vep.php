<?php 
/** 
	@page an_vep
	
	@todo test extended splice region plugin: https://github.com/Ensembl/VEP_plugins/blob/release/94/SpliceRegion.pm
	@todo test parameters: --gene_phenotype --ccds --biotype --canonical --pubmed
	@todo check bugs are fixed in new VEP version:
		- bug: chr2:48033890 CT>C not annotated with dbSNP identifier rs267608137
		- bug: chr17:41244936 G>A not annotated with 1000 Genomes AF
		- bug: chr5:138658646 A>C not annotated with dbSNP identifier, 1000 Genomes AF, genomAD AF
		- bug: chr16:89346631 G>C not annotated with dbSNP identifier, 1000 Genomes AF, genomAD AF
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_vep", "Variant annotation with Ensembl VEP.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFlag("all_transcripts", "Annotate all transcripts - if unset only GENCODE basic transcripts are annotated.");
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
extract($parser->parse($argv));

//annotate only field we really need to prevent bloatinf the VCF file 
$fields = array("Allele", "Consequence", "IMPACT", "SYMBOL", "Feature", "Feature_type", "EXON", "INTRON", "HGVSc", "HGVSp", "DOMAINS", "SIFT", "PolyPhen", "Existing_variation", "AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF", "EA_AF", "AA_AF"); 

$vep_path = dirname(get_path("vep"));
$vep_data_path = get_path("local_data")."/".basename(get_path("vep_data"))."/"; //the data is copied to the local data folder by 'data_setup' to speed up annotations (and prevent hanging annotation jobs)

$args = array();
$args[] = "-i $in --format vcf"; //input
$args[] = "-o $out --vcf --no_stats --force_overwrite"; //output
$args[] = "--species homo_sapiens --assembly {$build}"; //species
$args[] = "--fork {$threads}"; //speed (--buffer_size did not change run time when between 1000 and 20000)
$args[] = "--offline --cache --dir_cache {$vep_data_path}/ --fasta ".get_path("local_data")."/{$build}.fa"; //paths to data
$args[] = "--numbers --hgvs --domains"; //annotation options
$args[] = "--regulatory"; //regulatory features
$fields[] = "BIOTYPE"; 
$args[] = "--sift p --polyphen p"; //pathogenicity predictions
$args[] = "--af --af_gnomad --af_esp --failed 1"; //population frequencies
$args[] = "--plugin CADD,".get_path("data_folder")."/dbs/CADD/whole_genome_SNVs.tsv.gz,".get_path("data_folder")."/dbs/CADD/InDels.tsv.gz"; //CADD
$fields[] = "CADD_RAW";
$args[] = "--plugin REVEL,".get_path("data_folder")."/dbs/REVEL/revel_all_chromosomes.tsv.gz"; //REVEL
$fields[] = "REVEL";
$args[] = "--plugin FATHMM_MKL,".get_path("data_folder")."/dbs/fathmm-MKL/fathmm-MKL_Current.tab.gz"; //fathmm-MKL
$fields[] = "FATHMM_MKL_C";
$fields[] = "FATHMM_MKL_NC";
$args[] = "--plugin MaxEntScan,{$vep_path}/MaxEntScan/"; //MaxEntScan
$fields[] = "MaxEntScan_ref";
$fields[] = "MaxEntScan_alt";
$args[] = "--plugin GeneSplicer,{$vep_path}/GeneSplicer/bin/linux/genesplicer,{$vep_path}/GeneSplicer/human,context=50"; //GeneSplicer
$fields[] = "GeneSplicer";
$args[] = "--plugin dbscSNV,".get_path("data_folder")."/dbs/dbscSNV/dbscSNV1.1_GRCh37.txt.gz"; //dbscSNV
$fields[] = "ada_score";
$fields[] = "rf_score";
$args[] = "--custom ".get_path("data_folder")."/dbs/gnomAD/gnomAD_genome_r2.0.2.vcf.gz,gnomADg,vcf,exact,0,AF,Hom,Hemi"; //genomAD
$fields[] = "gnomADg_AF";
$fields[] = "gnomADg_Hom";
$fields[] = "gnomADg_Hemi";
$args[] = "--custom ".get_path("data_folder")."/dbs/RepeatMasker/RepeatMasker.bed.gz,REPEATMASKER,bed,overlap,0"; //RepeatMasker
$fields[] = "REPEATMASKER";
$args[] = "--custom ".get_path("data_folder")."/dbs/ClinVar/clinvar_20180805_converted.vcf.gz,CLINVAR,vcf,exact,0,DETAILS"; //ClinVar
$fields[] = "CLINVAR";
$fields[] = "CLINVAR_DETAILS";
$args[] = "--custom ".get_path("data_folder")."/dbs/phyloP/hg19.100way.phyloP100way.bw,PHYLOP,bigwig"; //phyloP
$fields[] = "PHYLOP";
$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed.gz"; //OMIM annotation (optional because of license)
if(file_exists($omim_file))
{
	$args[] = "--custom {$omim_file},OMIM,bed,overlap,0";
	$fields[] = "OMIM";
}
$hgmd_file = get_path("data_folder")."/dbs/HGMD/HGMD_PRO_2018_3_fixed.vcf.gz"; //HGMD annotation (optional because of license)
if(file_exists($hgmd_file))
{
	$args[] = "--custom {$hgmd_file},HGMD,vcf,exact,0,CLASS,MUT,GENE,PHEN";
	$fields[] = "HGMD";
	$fields[] = "HGMD_CLASS";
	$fields[] = "HGMD_MUT";
	$fields[] = "HGMD_GENE";
	$fields[] = "HGMD_PHEN";
}	
if (!$all_transcripts)
{
	$args[] = "--gencode_basic";
}
$args[] = "--fields ".implode(",", $fields);
putenv("PERL5LIB={$vep_path}/Bio/:{$vep_path}/cpan/lib/perl5/:".getenv("PERL5LIB"));
$parser->exec(get_path("vep"), implode(" ", $args), true);

?>
