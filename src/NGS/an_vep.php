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
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFlag("all_transcripts", "Annotate all transcripts - if unset only GENCODE basic transcripts are annotated.");
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addFlag("skip_ngsd", "Do not annotate the NGSD counts.");
$parser->addFlag("somatic", "Also annotate the NGSD somatic counts.");
$parser->addFlag("test", "Use limited constant NGSD VCF file from test folder for annotation.");
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



// custom annotation by VcfAnnotateFromVcf

// create config file
$config_file_path = $parser->tempFile(".config");
$config_file = fopen($config_file_path, 'w');


// add gnomAD annotation
fwrite($config_file, annotation_file_path("/dbs/gnomAD/gnomAD_genome_r2.1.1.vcf.gz")."\tgnomADg\tAF,Hom,Hemi\t\ttrue\n");


// add clinVar annotation
fwrite($config_file, annotation_file_path("/dbs/ClinVar/clinvar_20190503_converted.vcf.gz")."\tCLINVAR\tDETAILS\tID\n");


// add HGMD annotation
$hgmd_file = annotation_file_path("/dbs/HGMD/HGMD_PRO_2019_2_fixed.vcf.gz", true); //HGMD annotation (optional because of license)
if(file_exists($hgmd_file))
{
	fwrite($config_file, $hgmd_file."\tHGMD\tCLASS,MUT,GENE,PHEN\tID\n");
}


if (!$skip_ngsd)
{
	// add NGSD annotation

	// determine NGSD annotation file:
	if ($test)
	{
		$ngsd_som_file = repository_basedir()."/test/data/an_vep_NGSD_somatic.vcf.gz";
		$ngsd_file = repository_basedir()."/test/data/an_vep_NGSD_germline.vcf.gz";
	}
	else
	{
		$ngsd_som_file = annotation_file_path("/dbs/NGSD/NGSD_somatic.vcf.gz");
		$ngsd_file = annotation_file_path("/dbs/NGSD/NGSD_germline.vcf.gz");
	}
	

	// get disease group column name:
	$disease_group_column = "";
	if (db_is_enabled("NGSD"))
	{
		if ($somatic)
		{
			$db_conn = DB::getInstance("NGSD", false);
			$disease_groups = $db_conn->getEnum("sample", "disease_group");
			trigger_error("All somatic samples belongs to disease group Neoplasms", E_USER_NOTICE);
			$disease_group_column_idx = array_search("Neoplasms", $disease_groups);
			$disease_group_column = "GSC".sprintf('%02d', $disease_group_column_idx + 1);
		}
		else
		{
			// check file format
			if (ends_with($in, ".vcf"))
			{
				// parse vcf to get sample name
				$handle = fopen2($in, "r");
				$ps_name = "";
				while(!feof($handle))
				{
					
					$line = nl_trim(fgets($handle));
					if($line=="") continue;

					//header line
					if(starts_with($line, "#CHROM")) 
					{
						$parts = explode("\t", $line);
						$format_idx = array_search("FORMAT", $parts);
						$ps_name = $parts[($format_idx + 1)];
						trigger_error("Processed sample id: $ps_name ", E_USER_NOTICE);
						break;
					}
				}
				if ($ps_name != "")
				{
					$db_conn = DB::getInstance("NGSD", false);
					$disease_groups = $db_conn->getEnum("sample", "disease_group");
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
						trigger_error("No disease group found for sample $ps_name. NGSD count annotation for disease group will be missing in output file.",E_USER_WARNING);
					}
				}
				else
				{
					trigger_error("No sample name found in VCF file. NGSD count annotation for disease group will be missing in output file.",E_USER_WARNING);
				}
			}
			else
			{
				trigger_error("Can not extract disease group from gzipped VCF files. NGSD count annotation for disease group will be missing in output file.",E_USER_WARNING);
			}
		}
		
	}
	else
	{
		trigger_error("No connection to NGSD. NGSD count annotation for disease group will be missing in output file.",E_USER_WARNING);
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
		trigger_error("VCF file for NGSD germline annotation not found at '".$ngsd_file."'. NGSD annotation will be missing in output file.",E_USER_WARNING);
	}
	

	if ($somatic)
	{
		if (file_exists($ngsd_som_file))
		{
			fwrite($config_file, $ngsd_som_file."\tNGSD\tSOM_C,SOM_P\t\n");
		}
		else
		{
			trigger_error("VCF file for NGSD somatic annotation not found at '".$ngsd_som_file."'. NGSD annotation will be missing in output file.",E_USER_WARNING);
		}
	}
	
}


// close config file
fclose($config_file);


// execute VcfAnnotateFromVcf
$vafv_output = $parser->tempFile("_vafv.vcf");
$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromVcf", "-config_file ".$config_file_path." -in $vep_output -out $vafv_output", true);

if (!$skip_ngsd)
{
	// annotate genes
	$gene_file = annotation_file_path("/dbs/NGSD/NGSD_genes.bed");
	if (file_exists($gene_file))
	{
		$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromBed", "-bed ".$gene_file." -name NGSD_GENE_INFO -in $vafv_output -out $out", true);
	}
	else
	{
		trigger_error("BED file for NGSD gene annotation not found at '".$gene_file."'. NGSD annotation will be missing in output file.",E_USER_WARNING);
	}
}
else
{
	// use temp file as output
	$parser->moveFile($vafv_output, $out);
}


?>
