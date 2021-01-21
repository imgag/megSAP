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

//calculates MMSplice score for all variants of AF<=1
//annotates SpliceAI score for all variants inside NGSD + calculates score for private variants if less than threshold
function annotate_splice_predictions(&$vcf_annotate_output, $spliceai_threshold = 2000)
{
	global $build;
	global $threads;
	global $data_folder;
	global $out;
	global $parser;

	//annotate splice site variations with MMSplice for variants with AF <= 1 (1000genomes && gnomAD)
	$lowAF_variants = $parser->tempFile("_mmsplice_lowAF.vcf");
	$mmsplice_output = $parser->tempFile("_mmsplice.vcf");
	$gtf = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";
	$fasta = genome_fasta($build);
	$splice_env = get_path("Splicing", true);
	addMissingContigsToVcf($build, $vcf_annotate_output);
	$args = array();
	$args[] = "--vcf_in {$vcf_annotate_output}"; //input vcf
	$args[] = "--vcf_lowAF {$lowAF_variants}"; //input vcf storing only low allel frequency variants
	$args[] = "--vcf_out {$mmsplice_output}"; //output vcf
	$args[] = "--gtf {$gtf}"; //gtf annotation file
	$args[] = "--fasta {$fasta}"; //fasta reference file
	$args[] = "--threads {$threads}"; //fasta reference file
	putenv("PYTHONPATH");
	$parser->exec("OMP_NUM_THREADS={$threads} {$splice_env}/splice_env/bin/python3 ".repository_basedir()."/src/Tools/vcf_mmsplice_predictions.py", implode(" ", $args), true);

	// annotate SpliceAI for all NGSD variants
	$spliceai_file =  annotation_file_path("/dbs/SpliceAI/spliceai_scores.ngsd.13.12.20.vcf.gz");
	$spliceai_annotated_from_dbs = false;
	if (file_exists($spliceai_file))
	{
		$tmp = $parser->tempFile("_spai_preannotation.vcf");
		$parser->exec(get_path("ngs-bits") . "VcfAnnotateFromVcf", "-in $mmsplice_output -annotation_file $spliceai_file -info_ids SpliceAI -out $tmp" );
		$parser->moveFile($tmp, $mmsplice_output);
		$spliceai_annotated_from_dbs = true;
	}
	else
	{
		trigger_error("SpliceAI file for annotation not found at '".$spliceai_file."'. Private variants will be scored with SpliceAI if they do not exceed {$spliceai_threshold} variants.",E_USER_WARNING);
	}

	//store variants of AF<=1
	$private_var_dict = array();
	$low_af_file_mmsplice_h = fopen2($lowAF_variants, 'r');
	while(!feof($low_af_file_mmsplice_h))
	{
		$line = fgets($low_af_file_mmsplice_h);
		if (strlen(trim($line))==0) continue;
		//skip headers
		if(starts_with($line, "#"))
		{
			continue;
		}
		else
		{
			if(strlen(trim($line))==0) continue;
			$fields = explode("\t", $line);
			if(count($fields) < 5) continue;
			$var_ids = array($fields[0],$fields[1],$fields[3],$fields[4]);
			$private_var_dict[] = implode("", $var_ids); //string describing private variant
		}
	}
	fclose($low_af_file_mmsplice_h);

	//parse all lines that are not annotated with SpliceAI yet
	$wo_spliceai_annotation_file = $parser->tempFile("_wo_spliceai.vcf");
	exec2("grep '##fileformat' {$mmsplice_output} >> {$wo_spliceai_annotation_file}", true);
	exec2("grep '##contig' {$mmsplice_output} >> {$wo_spliceai_annotation_file}", true);
	exec2("grep '#CHROM' {$mmsplice_output} | cut -f1-8 -d'\t' >> {$wo_spliceai_annotation_file}", true);
	exec2("grep -v '#\|SpliceAI' {$mmsplice_output} | cut -f1-5 -d'\t' | sed 's/$/\t.\t.\t./' >> {$wo_spliceai_annotation_file}", false); //SpAI annotation might be empty

	//filter for regions of low AF
	$low_af_file_spliceai = $parser->tempFile("_spliceai_lowAF.vcf");
	$low_af_file_spliceai_h = fopen2($low_af_file_spliceai, 'w');
	$wo_spliceai_annotation_file_h = fopen2($wo_spliceai_annotation_file, 'r');
	while(!feof($wo_spliceai_annotation_file_h))
	{
		$line = fgets($wo_spliceai_annotation_file_h);
		//skip headers
		if(!starts_with($line, "chr"))
		{
			fwrite($low_af_file_spliceai_h, $line);
		}
		else
		{
			if(strlen(trim($line))==0) continue;
			$fields = explode("\t", $line);
			if(count($fields) < 5) continue;
			$var_ids = array($fields[0],$fields[1],$fields[3],$fields[4]);
			$needle = implode("", $var_ids);
			if(in_array($needle, $private_var_dict))
			{
				fwrite($low_af_file_spliceai_h, $line);
			}
		}
	}
	fclose($wo_spliceai_annotation_file_h);
	fclose($low_af_file_spliceai_h);

	//create bed file of SpliceAI variants that are scored 
	$new_spliceai_annotation = $parser->tempFile("_newSpliceAi_annotations.vcf");
	$spai_regions = $parser->tempFile("spai_scoring_regions.bed");
	$spai_build = strtolower($build);
	exec2("cut -f 2,4,5 -d'\t' {$splice_env}/splice_env/lib/python3.6/site-packages/spliceai/annotations/{$spai_build}.txt | sed 's/^/chr/' | sed '1d' > {$spai_regions}");
	$low_af_file_spliceai_filtered = $parser->tempFile("_private_spliceai.vcf");

	$parser->exec(get_path("ngs-bits")."/VcfFilter", "-reg ".$spai_regions." -in $low_af_file_spliceai -out $low_af_file_spliceai_filtered", true);
	list($private_variant_lines, $stderr)  = exec2("grep -v '##' $low_af_file_spliceai_filtered", false); //SpAI annotation might be empty
	$private_variant_count = count(array_filter($private_variant_lines));

	//calculate new SpliceAI score for private variants
	$spliceai_annotated = FALSE;
	if($private_variant_count <= $spliceai_threshold && $private_variant_count > 0)
	{
		$args = array();
		$args[] = "-I {$low_af_file_spliceai_filtered}";
		$args[] = "-O {$new_spliceai_annotation}";
		$args[] = "-R {$fasta}"; //output vcf
		$lower_build = strtolower($build);
		$args[] = "-A {$lower_build}"; //gtf annotation file
		putenv("PYTHONPATH");
		$parser->exec("OMP_NUM_THREADS={$threads} {$splice_env}/splice_env/bin/python3 {$splice_env}/splice_env/lib/python3.6/site-packages/spliceai", implode(" ", $args), true);
		$spliceai_annotated = TRUE;
	}
	else
	{
		trigger_error("SpliceAI annotation of private variants will be missing (More than ".$spliceai_threshold." variants).", E_USER_WARNING);
	}

	//annotate file with all splice predictions with calculated SpAI score for private variants
	if($spliceai_annotated)
	{
		//cut old header line, since new annotation produces additional line
		if($spliceai_annotated_from_dbs)
		{
			list($spai_header, $stderr)  = exec2("grep '##INFO=<ID=SpliceAI' {$mmsplice_output}");
			$spai_header = $spai_header[0]; //if grep produces no error there must be at least one match
			$spai_header = str_replace("\"", "\\\"", $spai_header); //SpAI header has doubles quotes which need to be escaped
			exec2("sed -i '/^##INFO=<ID=SpliceAI/d' ".$mmsplice_output);
		}

		//annotate new private variants
		$new_annotations_zipped = $parser->tempFile("_newSpliceAi_annotations.vcf.gz");
		$parser->exec("bgzip", "-c $new_spliceai_annotation > $new_annotations_zipped");
		$parser->exec("tabix", "-f -p vcf $new_annotations_zipped");
		$parser->exec(get_path("ngs-bits") . "VcfAnnotateFromVcf", "-in $mmsplice_output -annotation_file $new_annotations_zipped -info_ids SpliceAI -out $out" );

		//replace the header line with the origional one, station the annotation file of SpliceAI
		if($spliceai_annotated_from_dbs) exec2("sed -i 's/##INFO=<ID=SpliceAI.*/{$spai_header}/g' {$out}");
	}
	else
	{
		$parser->moveFile($mmsplice_output, $out);
	}
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
$args[] = "--numbers --hgvs --domains --transcript_version"; //annotation options
$args[] = "--regulatory"; //regulatory features
$fields[] = "BIOTYPE"; 
$args[] = "--sift b --polyphen b"; //pathogenicity predictions
$args[] = "--af --af_gnomad --failed 1"; //population frequencies
$args[] = "--plugin REVEL,".annotation_file_path("/dbs/REVEL/revel_all_chromosomes.tsv.gz"); //REVEL
$fields[] = "REVEL";
$args[] = "--plugin FATHMM_MKL,".annotation_file_path("/dbs/fathmm-MKL/fathmm-MKL_Current.tab.gz"); //fathmm-MKL
$fields[] = "FATHMM_MKL_C";
$fields[] = "FATHMM_MKL_NC";
$args[] = "--plugin MaxEntScan,{$vep_path}/MaxEntScan/"; //MaxEntScan
$fields[] = "MaxEntScan_ref";
$fields[] = "MaxEntScan_alt";
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

//annotate only fields we really need to prevent bloating the VCF file 
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


// custom annotation by VcfAnnotateFromVcf

// create config file
$config_file_path = $parser->tempFile(".config");
$config_file = fopen($config_file_path, 'w');


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

//add CADD score annotation
fwrite($config_file, annotation_file_path("/dbs/CADD/CADD_SNVs_1.6.vcf.gz")."\tCADD\tCADD=SNV\t\n");
fwrite($config_file, annotation_file_path("/dbs/CADD/CADD_InDels_1.6.vcf.gz")."\tCADD\tCADD=INDEL\t\n");

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
$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromVcf", "-config_file ".$config_file_path." -in $vep_output_refseq -out $vcf_annotate_output -threads $threads", true);

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
		$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromBed", "-bed ".$gene_file." -name NGSD_GENE_INFO -in $vcf_annotate_output -out $tmp", true);
		$parser->moveFile($tmp, $vcf_annotate_output);
	}
	else
	{
		trigger_error("BED file for NGSD gene annotation not found at '".$gene_file."'. NGSD annotation will be missing in output file.",E_USER_WARNING);
	}
}

//annotate MMSplice and SpliceAI predictions
annotate_splice_predictions($vcf_annotate_output);

//validate created VCF file
//check vcf file
if($check_lines >= 0)
{
	$parser->exec(get_path("ngs-bits")."VcfCheck", "-in $out -lines $check_lines -ref ".genome_fasta($build), true);
}

?>