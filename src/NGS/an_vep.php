<?php 
/** 
	@page an_vep
	
	@todo test Mastermind plugin: https://github.com/Ensembl/VEP_plugins/blob/release/98/Mastermind.pm
	@todo test FunMotifs plugin: https://github.com/Ensembl/VEP_plugins/blob/release/98/FunMotifs.pm
	@todo test parameters: --gene_phenotype --ccds --biotype --canonical
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_vep", "Variant annotation with Ensembl VEP.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
//optional
$parser->addString("ps_name", "Processed sample name (used to determine sample meta info from NGSD (e.g. disease group).", true, "");
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFlag("all_transcripts", "Annotate all transcripts - if unset only GENCODE basic transcripts are annotated.");
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addFlag("somatic", "Also annotate the NGSD somatic counts.");
$parser->addFlag("no_splice", "Skip splicing predictions of private variants (this can be very slow).");
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
$fields = array("Allele", "Consequence", "IMPACT", "SYMBOL", "HGNC_ID", "Feature", "Feature_type", "EXON", "INTRON", "HGVSc", "HGVSp", "DOMAINS", "SIFT", "PolyPhen", "Existing_variation", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF");

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
$args[] = "--numbers --hgvs --domains --transcript_version"; //annotation options
$args[] = "--regulatory"; //regulatory features
$fields[] = "BIOTYPE"; 
$args[] = "--sift b --polyphen b"; //pathogenicity predictions
$args[] = "--af_gnomad --failed 1"; //population frequencies
$args[] = "--plugin REVEL,".annotation_file_path("/dbs/REVEL/revel_grch38_all_chromosomes.tsv.gz"); //REVEL
$fields[] = "REVEL";
$args[] = "--plugin MaxEntScan,{$vep_path}/MaxEntScan/"; //MaxEntScan
$fields[] = "MaxEntScan_ref";
$fields[] = "MaxEntScan_alt";
$args[] = "--pubmed"; //add publications
$fields[] = "PUBMED";
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

// add phyloP annotation:
$vcf_output_bigwig = $parser->tempFile("_bigwig.vcf");
$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromBigWig", "-name PHYLOP -desc \"".annotation_file_path("/dbs/phyloP/hg38.phyloP100way.bw")." (ngs-bits/VcfAnnotateFromBigWig - mode max)\" -mode max -in {$vep_output} -out {$vcf_output_bigwig} -bw ".annotation_file_path("/dbs/phyloP/hg38.phyloP100way.bw")." -threads {$threads}", true);


// custom annotation by VcfAnnotateFromVcf

// create config file
$config_file_path = $parser->tempFile(".config");
$config_file = fopen2($config_file_path, 'w');


// add gnomAD annotation
fwrite($config_file, annotation_file_path("/dbs/gnomAD/gnomAD_genome_v3.1.2_GRCh38.vcf.gz")."\tgnomADg\tAC,AF,Hom,Hemi,Het,Wt\t\ttrue\n");
fwrite($config_file, annotation_file_path("/dbs/gnomAD/gnomAD_genome_v3.1.mito_GRCh38.vcf.gz")."\tgnomADm\tAF_hom\t\ttrue\n");

// add clinVar annotation
fwrite($config_file, annotation_file_path("/dbs/ClinVar/clinvar_20220328_converted_GRCh38.vcf.gz")."\tCLINVAR\tDETAILS\tID\n");

// add HGMD annotation
$hgmd_file = annotation_file_path("/dbs/HGMD/HGMD_PRO_2022_1_fixed.vcf.gz", true); //HGMD annotation (optional because of license)
if(file_exists($hgmd_file))
{
	fwrite($config_file, $hgmd_file."\tHGMD\tCLASS,MUT,GENE,PHEN\tID\n");
}

//add CADD score annotation
fwrite($config_file, annotation_file_path("/dbs/CADD/CADD_SNVs_1.6_GRCh38.vcf.gz")."\tCADD\tCADD=SNV\t\n");
fwrite($config_file, annotation_file_path("/dbs/CADD/CADD_InDels_1.6_GRCh38.vcf.gz")."\tCADD\tCADD=INDEL\t\n");

//precalculated SpliceAI scores
$spliceai_file = annotation_file_path("/dbs/SpliceAI/spliceai_scores_2022_05_24_GRCh38.vcf.gz");
if (file_exists($spliceai_file))
{
	fwrite($config_file, $spliceai_file."\t\tSpliceAI\t\n");
}
else
{
	trigger_error("SpliceAI annotation file with pre-calcualted scores not found at '".$spliceai_file."'.", E_USER_WARNING);
	$no_splice = true;
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
			fwrite($config_file, $ngsd_som_file."\tNGSD\tSOM_C,SOM_P,SOM_VICC,SOM_VICC_COMMENT\t\n");
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
$vcf_annotate_output = $parser->tempFile("_annotateFromVcf.vcf");
$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromVcf", "-config_file ".$config_file_path." -in {$vcf_output_bigwig} -out {$vcf_annotate_output} -threads {$threads}", true);

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
		$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromBed", "-bed ".$gene_file." -name NGSD_GENE_INFO -sep '&' -in {$vcf_annotate_output} -out {$tmp} -threads {$threads}", true);
		$parser->moveFile($tmp, $vcf_annotate_output);
	}
	else
	{
		trigger_error("BED file for NGSD gene annotation not found at '".$gene_file."'. NGSD annotation will be missing in output file.",E_USER_WARNING);
	}
}

//perform splicing predictions of private variants using SpliceAI (very slow)
if (!$no_splice)
{
	$tmp = $parser->tempFile("_spliceai_private.vcf");
	$parser->execTool("NGS/an_spliceai.php", "-in {$vcf_annotate_output} -out {$tmp} -threads {$threads} -build {$build}");
	$parser->moveFile($tmp, $vcf_annotate_output);
}

//annotate RepeatMasker
$tmp = $parser->tempFile("_repeatmasker.vcf");
$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromBed", "-bed ".annotation_file_path("/dbs/RepeatMasker/RepeatMasker_GRCh38.bed")." -name REPEATMASKER -sep '&' -in {$vcf_annotate_output} -out {$tmp} -threads {$threads}", true);
$parser->moveFile($tmp, $vcf_annotate_output);

//annotate OMIM (optional because of license)
$omim_file = annotation_file_path("/dbs/OMIM/omim.bed", true);
if(file_exists($omim_file))
{
	$tmp = $parser->tempFile("_omim.vcf");
	$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromBed", "-bed {$omim_file} -name OMIM -sep '&' -in {$vcf_annotate_output} -out {$tmp}  -threads {$threads}", true);
	$parser->moveFile($tmp, $vcf_annotate_output);
}

//mark variants in low-confidence regions
$low_conf_bed = repository_basedir()."/data/misc/low_conf_regions.bed";
$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $vcf_annotate_output -mark low_conf_region -inv -reg $low_conf_bed -out $out", true);

//validate created VCF file
//check vcf file
if($check_lines >= 0)
{
	$parser->exec(get_path("ngs-bits")."VcfCheck", "-in $out -lines $check_lines -ref ".genome_fasta($build), true);
}

?>