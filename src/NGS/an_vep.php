<?php 
/** 
	@page an_vep
	
	@todo test extended splice region plugin: https://github.com/Ensembl/VEP_plugins/blob/release/94/SpliceRegion.pm
	@todo test parameters: --gene_phenotype --ccds --biotype --canonical --pubmed
	@todo check annotation bugs are fixed in /mnt/users/ahsturm1/Sandbox/bugs/vep/annotations/
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

//get local/global data file path - depending on what is available
function annotation_file_path($rel_path, $is_optional=false)
{
	global $data_folder;
	global $local_data;
	
	//check if original exists
	$orig = $data_folder.$rel_path;
	if (!file_exists($orig))
	{
		if ($is_optional)
		{
			return $orig;
		}
		trigger_error("VEP annotation file '$orig' missing!", E_USER_ERROR);
	}
	
	//check if copy exists
	$copy = $local_data."/ensembl-vep-dbs/".basename($rel_path);
	if (!file_exists($copy))
	{
		trigger_error("VEP annotation file '$rel_path' not found in local data copy. Using (possibly slow) remote file '$orig'!", E_USER_NOTICE);
		return $orig;
	}
	
	return $copy;
}

//annotate only fields we really need to prevent bloating the VCF file 
$fields = array("Allele", "Consequence", "IMPACT", "SYMBOL", "HGNC_ID", "Feature", "Feature_type", "EXON", "INTRON", "HGVSc", "HGVSp", "DOMAINS", "SIFT", "PolyPhen", "Existing_variation", "AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF"); 

$vep_path = dirname(get_path("vep"));
$local_data = get_path("local_data");
$vep_data_path = "{$local_data}/".basename(get_path("vep_data"))."/"; //the data is copied to the local data folder by 'data_setup' to speed up annotations (and prevent hanging annotation jobs)
$data_folder = get_path("data_folder");

$args = array();
$args[] = "-i $in --format vcf"; //input
$args[] = "-o $out --vcf --no_stats --force_overwrite"; //output
$args[] = "--species homo_sapiens --assembly {$build}"; //species
$args[] = "--fork {$threads}"; //speed (--buffer_size did not change run time when between 1000 and 20000)
$args[] = "--offline --cache --dir_cache {$vep_data_path}/ --fasta ".genome_fasta($build); //paths to data
$args[] = "--numbers --hgvs --domains"; //annotation options
$args[] = "--regulatory"; //regulatory features
$fields[] = "BIOTYPE"; 
$args[] = "--sift b --polyphen b"; //pathogenicity predictions
$args[] = "--af --af_gnomad --af_esp --failed 1"; //population frequencies
$args[] = "--plugin CADD,".annotation_file_path("/dbs/CADD/whole_genome_SNVs.tsv.gz").",".annotation_file_path("/dbs/CADD/InDels.tsv.gz"); //CADD
$fields[] = "CADD_PHRED";
$args[] = "--plugin REVEL,".annotation_file_path("/dbs/REVEL/revel_all_chromosomes.tsv.gz"); //REVEL
$fields[] = "REVEL";
$args[] = "--plugin FATHMM_MKL,".annotation_file_path("/dbs/fathmm-MKL/fathmm-MKL_Current.tab.gz"); //fathmm-MKL
$fields[] = "FATHMM_MKL_C";
$fields[] = "FATHMM_MKL_NC";
$args[] = "--plugin MaxEntScan,{$vep_path}/MaxEntScan/"; //MaxEntScan
$fields[] = "MaxEntScan_ref";
$fields[] = "MaxEntScan_alt";
$args[] = "--plugin GeneSplicer,{$vep_path}/GeneSplicer/sources/genesplicer,{$local_data}/GeneSplicer/,tmpdir=".sys_get_temp_dir(); //GeneSplicer
$fields[] = "GeneSplicer";
$args[] = "--plugin dbscSNV,".annotation_file_path("/dbs/dbscSNV/dbscSNV1.1_GRCh37.txt.gz"); //dbscSNV
$fields[] = "ada_score";
$fields[] = "rf_score";
$args[] = "--custom ".annotation_file_path("/dbs/gnomAD/gnomAD_genome_r2.1.1.vcf.gz").",gnomADg,vcf,exact,0,AF,Hom,Hemi"; //genomAD
$fields[] = "gnomADg_AF";
$fields[] = "gnomADg_Hom";
$fields[] = "gnomADg_Hemi";
$args[] = "--custom ".annotation_file_path("/dbs/RepeatMasker/RepeatMasker.bed.gz").",REPEATMASKER,bed,overlap,0"; //RepeatMasker
$fields[] = "REPEATMASKER";
$args[] = "--custom ".annotation_file_path("/dbs/ClinVar/clinvar_20190503_converted.vcf.gz").",CLINVAR,vcf,exact,0,DETAILS"; //ClinVar
$fields[] = "CLINVAR";
$fields[] = "CLINVAR_DETAILS";
$args[] = "--custom ".annotation_file_path("/dbs/phyloP/hg19.100way.phyloP100way.bw").",PHYLOP,bigwig"; //phyloP
$fields[] = "PHYLOP";
$omim_file = annotation_file_path("/dbs/OMIM/omim.bed.gz", true); //OMIM annotation (optional because of license)
if(file_exists($omim_file))
{
	$args[] = "--custom {$omim_file},OMIM,bed,overlap,0";
	$fields[] = "OMIM";
}
$hgmd_file = annotation_file_path("/dbs/HGMD/HGMD_PRO_2019_2_fixed.vcf.gz", true); //HGMD annotation (optional because of license)
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

//print VEP warnings
$warn_file = $out."_warnings.txt";
if (file_exists($warn_file))
{
	$file = file($warn_file);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		print $line."\n";
	}
}

?>
