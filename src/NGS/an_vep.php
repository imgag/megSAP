<?php
/**
	@page an_vep

	@todo test extended splice region plugin: https://github.com/Ensembl/VEP_plugins/blob/release/94/SpliceRegion.pm
	@todo test Mastermind plugin: https://github.com/Ensembl/VEP_plugins/blob/release/98/Mastermind.pm
	@todo test FunMotifs plugin: https://github.com/Ensembl/VEP_plugins/blob/release/98/FunMotifs.pm
	@todo test parameters: --gene_phenotype --ccds --biotype --canonical --pubmed
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_vep", "Variant annotation with Ensembl VEP.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
//optional
$parser->addString("ps_name", "Processed sample name (used to determine sample meta info from NGSD (e.g. disease group).", true, "");
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFlag("all_transcripts", "Annotate all transcripts - if unset only GENCODE basic transcripts are annotated.");
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addFlag("somatic", "Also annotate the NGSD somatic counts.");
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

// generate temp file for vep output
$vep_output = $parser->tempFile("_vep.vcf");

//annotate only fields we really need to prevent bloating the VCF file
$fields = array("Allele", "Consequence", "IMPACT", "SYMBOL", "HGNC_ID", "Feature", "Feature_type", "EXON", "INTRON", "HGVSc", "HGVSp", "DOMAINS", "SIFT", "PolyPhen", "Existing_variation", "AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF");

$vep_path = dirname(get_path("vep"));
$local_data = get_path("local_data");
$vep_data_path = "{$local_data}/".basename(get_path("vep_data"))."/"; //the data is copied to the local data folder by 'data_setup' to speed up annotations (and prevent hanging annotation jobs)
$data_folder = get_path("data_folder");

$args = array();
$args[] = "-i $in --format vcf"; //input
$args[] = "-o $vep_output --vcf --no_stats --force_overwrite"; //output
$args[] = "--species homo_sapiens --assembly {$build}"; //species
$args[] = "--fork {$threads}"; //speed (--buffer_size did not change run time when between 1000 and 20000)
$args[] = "--offline --cache --dir_cache {$vep_data_path}/ --fasta ".genome_fasta($build); //paths to data
$args[] = "--numbers --hgvs --domains"; //annotation options
$args[] = "--regulatory"; //regulatory features
$fields[] = "BIOTYPE";
$args[] = "--sift b --polyphen b"; //pathogenicity predictions
$args[] = "--af --af_gnomad --failed 1"; //population frequencies
$args[] = "--plugin CADD,".annotation_file_path("/dbs/CADD/CADD_SNVs_1.6.tsv.gz").",".annotation_file_path("/dbs/CADD/CADD_InDels_1.6.tsv.gz"); //CADD
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
$fields[] = "GeneSplicer"; //TODO reactivate GeneSplicer when runtime problems are fixed in VEP 100 (https://github.com/Ensembl/ensembl-vep/issues/641)
$args[] = "--plugin dbscSNV,".annotation_file_path("/dbs/dbscSNV/dbscSNV1.1_GRCh37.txt.gz"); //dbscSNV
$fields[] = "ada_score";
$fields[] = "rf_score";;
$args[] = "--custom ".annotation_file_path("/dbs/RepeatMasker/RepeatMasker.bed.gz").",REPEATMASKER,bed,overlap,0"; //RepeatMasker
$fields[] = "REPEATMASKER";
$args[] = "--custom ".annotation_file_path("/dbs/phyloP/hg19.100way.phyloP100way.bw").",PHYLOP,bigwig"; //phyloP
$fields[] = "PHYLOP";

$omim_file = annotation_file_path("/dbs/OMIM/omim.bed.gz", true); //OMIM annotation (optional because of license)

if(file_exists($omim_file))
{
	$args[] = "--custom {$omim_file},OMIM,bed,overlap,0";
	$fields[] = "OMIM";
}

if (!$all_transcripts)
{
	$args[] = "--gencode_basic";
}

$args[] = "--fields ".implode(",", $fields);
putenv("PERL5LIB={$vep_path}/Bio/:{$vep_path}/cpan/lib/perl5/:".getenv("PERL5LIB"));
$parser->exec(get_path("vep"), implode(" ", $args), true);

//print VEP warnings
$warn_file = $vep_output."_warnings.txt";
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

// second VEP run to annotate RefSeq transcripts
// generate temp file for vep output
$vep_output_refseq = $parser->tempFile("_vep.vcf");

// annotate only fields we really need to prevent bloating the VCF file
$fields = array("Allele", "Consequence", "IMPACT", "SYMBOL", "HGNC_ID", "Feature", "Feature_type", "EXON", "INTRON", "HGVSc", "HGVSp", "DOMAINS");

$args = array();
$args[] = "-i $vep_output --format vcf"; //input
$args[] = "-o $vep_output_refseq --vcf --no_stats --force_overwrite"; //output
$args[] = "--species homo_sapiens --assembly {$build}"; //species
$args[] = "--refseq"; //use RefSeq annotation instead of ensemble
$args[] = "--vcf_info_field CSQ_refseq"; //store annotation in custom info field to keep RefSeq and ensemble separate
$args[] = "--fork {$threads}"; //speed (--buffer_size did not change run time when between 1000 and 20000)
$args[] = "--offline --cache --dir_cache {$vep_data_path}/ --fasta ".genome_fasta($build); //paths to data
$args[] = "--numbers --hgvs --domains"; //annotation options


$args[] = "--fields ".implode(",", $fields);
putenv("PERL5LIB={$vep_path}/Bio/:{$vep_path}/cpan/lib/perl5/:".getenv("PERL5LIB"));
$parser->exec(get_path("vep"), implode(" ", $args), true);

//print VEP warnings
$warn_file = $vep_output_refseq."_warnings.txt";
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

// TODO: prefilter VCF to only pass coding variants to the AIdiva subroutine
$coding_filtered = $parser->tempFile("_coding.vcf");
$parser->exec("python3 ".get_path("aidiva")."aidiva/helper_modules/filter_vcf.py", "--in_file $vep_output --out_file $coding_filtered", true);


$family_file = ""; // handle sample as single sample
$aidiva_config = get_path("aidiva")."/data/AIdiva_configuration_smallTestFile_annotated.yaml";
$ref_genome = genome_fasta($build);

$temp_results = $parser->tempFolder("aidiva_workdir");
$args = array();
$args[] = "-vcf {$coding_filtered}";
$args[] = "-outdir {$temp_results}";
if ($family_file != "")
{
	$args[] = "-family {$family_file}";
}
if ($ps_name != "")
{
	$args[] = "-ps_name {$ps_name}";
}
$args[] = "-genome_file {$ref_genome}";
$args[] = "-config {$aidiva_config}";
$args[] = "-threads {$threads}";
$parser->execTool("NGS/sp_aidiva.php", implode(" ", $args));

$aidiva_result_file = $temp_results."/"."aidiva"."_result_sorted.vcf.gz";


// custom annotation by VcfAnnotateFromVcf

// create config file
$config_file_path = $parser->tempFile(".config");
$config_file = fopen($config_file_path, 'w');

// add AIdiva annotation
fwrite($config_file, $aidiva_result_file."\t\tAIDIVA\t\ttrue\n");


// add gnomAD annotation
fwrite($config_file, annotation_file_path("/dbs/gnomAD/gnomAD_genome_r2.1.1.vcf.gz")."\tgnomADg\tAF,Hom,Hemi\t\ttrue\n");


// add clinVar annotation
fwrite($config_file, annotation_file_path("/dbs/ClinVar/clinvar_20200506_converted.vcf.gz")."\tCLINVAR\tDETAILS\tID\n");


// add HGMD annotation
$hgmd_file = annotation_file_path("/dbs/HGMD/HGMD_PRO_2020_1_fixed.vcf.gz", true); //HGMD annotation (optional because of license)
if(file_exists($hgmd_file))
{
	fwrite($config_file, $hgmd_file."\tHGMD\tCLASS,MUT,GENE,PHEN\tID\n");
}

// check if NGSD export file is available:
$skip_ngsd = false;
$ngsd_file = $data_folder."/dbs/NGSD/NGSD_germline.vcf.gz";
$ngsd_som_file = $data_folder."/dbs/NGSD/NGSD_somatic.vcf.gz";
if (!file_exists($ngsd_file))
{
	trigger_error("VCF file for NGSD germline annotation not found at '".$ngsd_file."'. NGSD annotation will be missing in output file.",E_USER_WARNING);
	$skip_ngsd = true;
}
if ($somatic && !file_exists($ngsd_som_file))
{
	trigger_error("VCF file for NGSD somatic annotation not found at '".$ngsd_som_file."'. NGSD annotation will be missing in output file.",E_USER_WARNING);
	$skip_ngsd = true;
}

if (!$skip_ngsd)
{
	// add NGSD annotation
	if ($test)
	{
		$ngsd_som_file = repository_basedir()."/test/data/an_vep_NGSD_somatic.vcf.gz";
		$ngsd_file = repository_basedir()."/test/data/an_vep_NGSD_germline.vcf.gz";
	}

	// get disease group column name
	$disease_group_column = "";
	if (db_is_enabled("NGSD"))
	{
		$db_conn = DB::getInstance("NGSD", false);
		$disease_groups = $db_conn->getEnum("sample", "disease_group");
		if ($somatic)
		{
			trigger_error("All somatic samples belongs to disease group Neoplasms", E_USER_NOTICE);
			$disease_group_column_idx = array_search("Neoplasms", $disease_groups);
			$disease_group_column = "GSC".sprintf('%02d', $disease_group_column_idx + 1);
		}
		else //germline
		{
			if ($ps_name!="")
			{
				$details = get_processed_sample_info($db_conn, $ps_name, false);
				$disease_group = $details['disease_group'];
				$disease_group_column = "";
				if ($disease_group != "")
				{
					trigger_error("Sample '$ps_name' belongs to disease group $disease_group.", E_USER_NOTICE);
					$disease_group_column_idx = array_search($disease_group, $disease_groups);
					$disease_group_column = "GSC".sprintf('%02d', $disease_group_column_idx + 1);
				}
				else
				{
					trigger_error("NGSD count annotation for disease group will be missing in output file (no disease group set in NGSD for sample '$ps_name').", E_USER_WARNING);
				}
			}
			else
			{
				trigger_error("NGSD count annotation for disease group will be missing in output file ('ps_name' not given).", E_USER_WARNING);
			}
		}
	}
	else
	{
		trigger_error("NGSD count annotation for disease group will be missing in output file (NGSD is disabled).", E_USER_WARNING);
	}

	$ngsd_columns = ["COUNTS"];
	if ($disease_group_column != "")
	{
		$ngsd_columns[] = $disease_group_column."=GROUP";
	}
	array_push($ngsd_columns, "HAF", "CLAS", "CLAS_COM", "COM");

	if (file_exists($ngsd_file))
	{
		fwrite($config_file, $ngsd_file."\tNGSD\t".implode(",", $ngsd_columns)."\t\n");
	}
	else
	{
		trigger_error("VCF file for NGSD germline annotation not found at '".$ngsd_file."'!",E_USER_ERROR);
	}


	if ($somatic)
	{
		if (file_exists($ngsd_som_file))
		{
			fwrite($config_file, $ngsd_som_file."\tNGSD\tSOM_C,SOM_P\t\n");
		}
		else
		{
			trigger_error("VCF file for NGSD somatic annotation not found at '".$ngsd_som_file."'!",E_USER_ERROR);
		}

		// store file date of NGSD files to detect file changes during annotation
		$ngsd_som_file_mtime = filemtime($ngsd_som_file);
		if ($ngsd_som_file_mtime == false)
		{
			trigger_error("Cannot get modification date of '".$ngsd_som_file."'!",E_USER_ERROR);
		}
	}

	// store file date of NGSD files to detect file changes during annotation
	$ngsd_file_mtime = filemtime($ngsd_file);
	if ($ngsd_file_mtime == false)
	{
		trigger_error("Cannot get modification date of '".$ngsd_file."'!",E_USER_ERROR);
	}
}


// close config file
fclose($config_file);

// execute VcfAnnotateFromVcf
$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromVcf", "-config_file ".$config_file_path." -in $vep_output_refseq -out $out -threads $threads", true);

if (!$skip_ngsd)
{
	// check if files have changed during annotation:
	if (!($ngsd_file_mtime == filemtime($ngsd_file)))
	{
		trigger_error("Annotation file '".$ngsd_file."' has changed during annotation!",E_USER_ERROR);
	}
	if ($somatic)
	{
		if (!($ngsd_som_file_mtime == filemtime($ngsd_som_file)))
		{
			trigger_error("Annotation file '".$ngsd_som_file."' has changed during annotation!",E_USER_ERROR);
		}
	}

	// annotate genes
	$gene_file = $data_folder."/dbs/NGSD/NGSD_genes.bed";
	if (file_exists($gene_file))
	{

		$tmp = $parser->tempFile(".vcf");
		$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromBed", "-bed ".$gene_file." -name NGSD_GENE_INFO -in $out -out $tmp", true);
		$parser->moveFile($tmp, $out);
	}
	else
	{
		trigger_error("BED file for NGSD gene annotation not found at '".$gene_file."'. NGSD annotation will be missing in output file.",E_USER_WARNING);
	}
}



//validate created VCF file

//check vcf file
if($check_lines >= 0)
{
	$parser->exec(get_path("ngs-bits")."VcfCheck", "-in $out -lines $check_lines -ref ".genome_fasta($build), true);
}

$parser->exec("bgzip", "$out", true);

?>
