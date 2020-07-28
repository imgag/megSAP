<?php 
/** 
	@page an_aidiva_vep
	
	@todo remove all unnecessary annotations
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_aidiva_vep", "Variant annotation with Ensembl VEP.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFlag("all_transcripts", "Annotate all transcripts - if unset only GENCODE basic transcripts are annotated.");
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addFlag("somatic", "Also annotate the NGSD somatic counts.");
$parser->addFlag("basic", "Skips major part of the annotation.");
$parser->addFlag("test", "Use limited constant NGSD VCF file from test folder for annotation.");
$parser->addInt("check_lines", "Number of VCF lines that will be validated in the output file. (If set to 0 all lines will be checked, if set to -1 the validation will be skipped.)", true, 1000);
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
$args[] = "--max_af --af --af_gnomad --failed 1"; //population frequencies
$fields[] = "MAX_AF";

if(!$basic)
{
$args[] = "--sift s --polyphen s"; //pathogenicity predictions
$args[] = "--plugin CADD,".annotation_file_path("/dbs/CADD/CADD_SNVs_1.6.tsv.gz").",".annotation_file_path("/dbs/CADD/CADD_InDels_1.6.tsv.gz"); //CADD
$fields[] = "CADD_PHRED";
$args[] = "--plugin REVEL,".annotation_file_path("/dbs/REVEL/revel_all_chromosomes.tsv.gz"); //REVEL
$fields[] = "REVEL";
$args[] = "--custom ".annotation_file_path("/dbs/UCSC/hg19.simpleRepeats.bedGraph.gz").",simpleRepeat,bed,overlap,0";
$fields[] = "simpleRepeat";
$args[] = "--custom ".annotation_file_path("/dbs/UCSC/hg19.segmentDuplication.bedGraph.gz").",segmentDuplication,bed,overlap,0";
$fields[] = "segmentDuplication";
$args[] = "--custom ".annotation_file_path("/dbs/ABB/hg19.abb_score.bedGraph.gz").",ABB_SCORE,bed,exact,0"; 
$fields[] = "ABB_SCORE";
$args[] = "--custom ".annotation_file_path("/dbs/Eigen/Eigen_phred_hg19_chrom1-22.vcf.gz").",eigen,vcf,exact,0,EIGEN_PHRED"; 
$fields[] = "custom_EIGEN_PHRED";
$args[] = "--custom ".annotation_file_path("/dbs/Condel/hg19.fannsdb_onlyCONDEL.vcf.gz").",fannsdb,vcf,exact,0,CONDEL"; 
$fields[] = "fannsdb_CONDEL";
$args[] = "--custom ".annotation_file_path("/dbs/fathmm-XF/hg19.fathmm_XF_onlyCoding.vcf.gz").",fathmm_xf,vcf,exact,0,FATHMM_XF"; 
$fields[] = "custom_FATHMM_XF";
$args[] = "--custom ".annotation_file_path("/dbs/MutationAssessor/MA_scores_rel3_hg19_full.vcf.gz").",mutationassessor,vcf,exact,0,MutationAssessor";
$fields[] = "custom_MutationAssessor";
$args[] = "--custom ".annotation_file_path("/dbs/phastCons/hg19.phastCons46way_mammal.bw").",phastCons46mammal,bigwig"; //
$fields[] = "phastCons46mammal";
$args[] = "--custom ".annotation_file_path("/dbs/phastCons/hg19.phastCons46way_primate.bw").",phastCons46primate,bigwig"; //
$fields[] = "phastCons46primate";
$args[] = "--custom ".annotation_file_path("/dbs/phastCons/hg19.phastCons46way_vertebrate.bw").",phastCons46vertebrate,bigwig"; //
$fields[] = "phastCons46vertebrate";
$args[] = "--custom ".annotation_file_path("/dbs/phyloP/hg19.phyloP46way_mammal.bw").",phyloP46mammal,bigwig"; //
$fields[] = "phyloP46mammal";
$args[] = "--custom ".annotation_file_path("/dbs/phyloP/hg19.phyloP46way_primate.bw").",phyloP46primate,bigwig"; //
$fields[] = "phyloP46primate";
$args[] = "--custom ".annotation_file_path("/dbs/phyloP/hg19.phyloP46way_vertebrate.bw").",phyloP46vertebrate,bigwig"; //
$fields[] = "phyloP46vertebrate";
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

//validate created VCF file

//check vcf file
if($check_lines >= 0)
{
	$parser->exec(get_path("ngs-bits")."VcfCheck", "-in $out -lines $check_lines -ref ".genome_fasta($build), true);
}

?>
