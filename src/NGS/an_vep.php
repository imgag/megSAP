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
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFlag("all_transcripts", "Annotate all transcripts - if unset only GENCODE basic transcripts are annotated.");
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addFlag("somatic", "Also annotate the NGSD somatic counts.");
$parser->addFlag("no_splice", "Skip splicing predictions (SpliceAI).");
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

function vcf_variant_count($vcf)
{
	$count = 0;
	$h = fopen2($vcf, 'r');
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || $line[0]=='#') continue;
		
		++$count;
	}
	fclose($h);
	
	return $count;
}

function write_lowAF_variants($vcf_in, $max_af, $anno_key, $vcf_out)
{
	//init
	$csq_tg_idx = -1;
	$csq_gnomad_idx = -1;
	$c_written = 0;
	
	$h = fopen2($vcf_in, 'r');
	$h_o = fopen2($vcf_out, 'w');
	while(!feof($h))
	{
		$line = nl_trim(fgets($h));
		if ($line=="") continue;
		
		//header
		if(starts_with($line, "#"))
		{
			fwrite($h_o, $line."\n");
			
			//extract VEP entry indices
			if(starts_with($line, '##INFO=<ID=CSQ,'))
			{
                preg_match('/.*Description=\"(.*)\".*/', $line, $match);
				if(count($match) < 2) continue;
                $entries = explode("|", $match[1]);
				for($i=0; $i<count($entries); ++$i)
				{
					$entry = trim($entries[$i]);
					if($entry=='AF')
					{
                        $csq_tg_idx = $i;
					}
					else if($entry=='gnomAD_AF')
					{
						$csq_gnomad_idx = $i;
					}
				}
			}	
		}
		else //content
		{
			$parts = explode("\t", $line);
			if(count($parts) < 8) trigger_error("Invalid VCF content line in {$vcf_in}: Missing INFO column.", E_USER_ERROR);
			
			//skip already annotated variants
			if (contains($line, $anno_key)) continue;
			
			//skip variant with too high AF
			$af = 0.0;
            $info = explode(';', $parts[7]);
			foreach($info as $info_entry)
			{
				if(starts_with($info_entry, 'CSQ='))
				{
					$info_entry = explode('=', $info_entry, 2);
					$csq_annotations = explode('|', $info_entry[1]);

					//1000G AF
					$value = $csq_annotations[$csq_tg_idx];
					if (is_numeric($value))
					{
						$af = max($af, $value);
					}
					//genomAD AF (exome)
					$value = $csq_annotations[$csq_gnomad_idx];
					if (is_numeric($value))
					{
						$af = max($af, $value);
					}
				}
				//gnomAD AF (genome)
				else if(starts_with($info_entry, 'gnomADg_AF='))
				{
					$info_entry = explode('=', $info_entry, 2);
					$value = $info_entry[1];
					if (is_numeric($value))
					{
						$af = max($af, $value);
					}
				}
				//gnomAD AF (mito)
				else if(starts_with($info_entry, 'gnomADm_AF_hom='))
				{
					$info_entry = explode('=', $info_entry, 2);
					$value = $info_entry[1];
					if (is_numeric($value))
					{
						$af = max($af, $value);
					}
				}
			}
			if($af>$max_af) continue;
			
			//write
			fwrite($h_o, $line."\n");
			++$c_written;
		}
	}
	fclose($h);
	fclose($h_o);
	
	return $c_written;
}

//checks if any contig line is given, if not adds all contig lines from reference genome
function add_missing_contigs($vcf)
{
	global $parser;
	global $build;
	
	//check if contig lines are contained
	$contains_contig_lines = false;
	$line_below_reference_info = -1;

	$file = fopen2($vcf, 'r');
	$line_nr = 0;
	while(!feof($file))
	{
		++$line_nr;
		$line = trim(fgets($file));

		if(starts_with($line, "##reference"))
		{
			$line_below_reference_info = $line_nr;
		}
		else if (starts_with($line, "##"))
		{
			if(starts_with($line, "##contig"))
			{
				$contains_contig_lines = true;
			}
		}
		else break;
	}
	fclose($file);
		
	//add contig lines if necessary
	if($contains_contig_lines)
	{
		$parser->log("Contig lines are present > nothing to add!");
	}
	else
	{
		$parser->log("Contig lines are missing > adding them from genome FASTA!");
		$new_contigs = array();
		$new_file_lines = file($vcf);
		$build_path = genome_fasta($build);
		list($chr_lines) = exec2("grep chr {$build_path}");
		foreach($chr_lines as $line)
		{
			$chr = "";
			$len = "";
			$parts = explode(" ", $line);
			if(sizeof($parts) >= 3)
			{
				//get chromosome
				$chr = $parts[0];
				$chr = ltrim($chr, '>');

				//get length
				preg_match('/.*:(\d+):.*$/', $parts[2], $matches);
				if(sizeof($matches)>=2)
				{
					$len = $matches[1];
				}
				else
				{
					// try to parse chr header in GRCH38 format
					foreach ($parts as $part) 
					{
						if(starts_with($part, "LN:"))
						{
							$len = explode(":", $part)[1];
						}
					}
					if (!isset($len)) trigger_error("No length information for chromosome ".$chr." found in line(".$line.") for reference genome(".$build.")", E_USER_WARNING);
				}

				//add information to contig array
				if($chr && $len)
				{
					$new_contigs[] = "##contig=<ID={$chr}, length={$len}>";
				}
			}
			else
			{
				trigger_error("Contig information could not be parsed for line(".$line.") from reference genome(".$build.").", E_USER_WARNING);
			}
		}

		if(empty($new_contigs))
		{
			trigger_error("No new contig lines were written for ".$vcf, E_USER_WARNING);
		}
		else
		{
			//write contigs in second line if no reference genome line is given
			if($line_below_reference_info==-1)
			{
				$line_below_reference_info = 1;
			}
			
			array_splice($new_file_lines, $line_below_reference_info, 0, $new_contigs);  
			file_put_contents($vcf, implode("\n", $new_file_lines));
		}
	}
}

function annotate_spliceai_score($splicing_output, $threshold = 15000)
{
	global $build;
	global $threads;
	global $data_folder;
	global $parser;

	//annotate precalculated SpliceAI scores
	$spliceai_file = annotation_file_path("/dbs/SpliceAI/spliceai_scores_2022_02_09_GRCh38.vcf.gz");
	$spliceai_annotated_from_dbs = false;
	if (file_exists($spliceai_file))
	{
		$tmp = $parser->tempFile("_spai_preannotation.vcf");
		$parser->exec(get_path("ngs-bits") . "VcfAnnotateFromVcf", "-in {$splicing_output} -annotation_file {$spliceai_file} -info_ids SpliceAI -out {$tmp}  -threads {$threads}");
		$parser->moveFile($tmp, $splicing_output);
		$spliceai_annotated_from_dbs = true;
	}
	else
	{
		trigger_error("SpliceAI file for annotation not found at '".$spliceai_file."'.", E_USER_WARNING);
	}

	//filter for private variants (low AF and not annotated in the step before)
	$low_af_file_spliceai = $parser->tempFile("_spliceai_lowAF.vcf");
	write_lowAF_variants($splicing_output, 0.01, "SpliceAI=", $low_af_file_spliceai);

	//filter based on SpliceAI trascript regions
	$spai_regions = $parser->tempFile("spai_scoring_regions.bed");
	$splice_env = get_path("Splicing", true);
	exec2("cut -f 2,4,5 -d'\t' {$splice_env}/splice_env/lib/python3.6/site-packages/spliceai/annotations/".strtolower($build).".txt | sed 's/^/chr/' | sed '1d' > {$spai_regions}");
	$low_af_file_spliceai_filtered = $parser->tempFile("_private_spliceai.vcf");
	$parser->exec(get_path("ngs-bits")."/VcfFilter", "-reg {$spai_regions} -in {$low_af_file_spliceai} -out {$low_af_file_spliceai_filtered}", true);
	$private_variant_count = vcf_variant_count($low_af_file_spliceai_filtered);
	
	//skip if too many private variants
	if($private_variant_count > $threshold)
	{
		trigger_error("SpliceAI annotation of private variants skipped: {$private_variant_count} variants is too many (threshold is {$threshold})!", E_USER_WARNING);
		return;
	}
	
	//calculate new SpliceAI score for private variants
	$parser->log("Scoring {$private_variant_count} variants with SpliceAI");
	$args = array();
	
	$low_af_file_spliceai_gzipped =  $parser->tempFile("_newSpliceAi_annotations.vcf.gz");
	$parser->exec(get_path("ngs-bits")."VcfSort", "-in {$low_af_file_spliceai_filtered} -out {$low_af_file_spliceai}");
	$parser->exec("bgzip", "-c $low_af_file_spliceai > $low_af_file_spliceai_gzipped");
	$parser->exec("tabix", "-f -p vcf $low_af_file_spliceai_gzipped");
	
	$new_spliceai_annotation = $parser->tempFile("_newSpliceAi_annotations.vcf");
	$args[] = "-I {$low_af_file_spliceai_gzipped}";
	$args[] = "-O {$new_spliceai_annotation}";
	$fasta = genome_fasta($build);
	$args[] = "-R {$fasta}"; //output vcf
	$lower_build = strtolower($build);
	$args[] = "-A {$lower_build}"; //gtf annotation file
	putenv("PYTHONPATH");
	$parser->exec("OMP_NUM_THREADS={$threads} {$splice_env}/splice_env/bin/python3 {$splice_env}/splice_env/lib/python3.6/site-packages/spliceai", implode(" ", $args), true);
	
	//annotate file with all splice predictions with calculated SpliceAI score for private variants
	$spai_anno_file = $parser->tempFile("_spliceaiPrivate.vcf");
	exec2("grep '#' {$new_spliceai_annotation} >> {$spai_anno_file}", true);
	exec2("grep -v '#' {$new_spliceai_annotation} | grep 'SpliceAI' >> {$spai_anno_file}", false);
	$private_variant_count = vcf_variant_count($spai_anno_file);
	if($private_variant_count > 0)
	{
		//cut old header line, since new annotation produces additional line
		if($spliceai_annotated_from_dbs)
		{
			list($spai_header, $stderr)  = exec2("grep '##INFO=<ID=SpliceAI' {$splicing_output}");
			$spai_header = $spai_header[0]; //if annotated from dbs, there must be one match
			$spai_header = str_replace("\"", "\\\"", $spai_header); //SpliceAI header has doubles quotes which need to be escaped
			exec2("sed -i '/^##INFO=<ID=SpliceAI/d' ".$splicing_output);
		}

		//annotate new private variants
		$new_annotations_zipped = $parser->tempFile("_newSpliceAi_annotations.vcf.gz");
		$parser->exec(get_path("ngs-bits")."VcfSort", "-in {$new_spliceai_annotation} -out {$new_spliceai_annotation}");
		$parser->exec("bgzip", "-c $new_spliceai_annotation > $new_annotations_zipped");
		$parser->exec("tabix", "-f -p vcf $new_annotations_zipped");
		$tmp = $parser->tempFile("_final_annotation.vcf.gz");
		$parser->exec(get_path("ngs-bits") . "VcfAnnotateFromVcf", "-in {$splicing_output} -annotation_file {$new_annotations_zipped} -info_ids SpliceAI -out {$tmp}  -threads {$threads}");
		$parser->moveFile($tmp, $splicing_output);

		//replace the header line with the origional one, station the annotation file of SpliceAI
		if($spliceai_annotated_from_dbs) exec2("sed -i 's/##INFO=<ID=SpliceAI.*/{$spai_header}/g' {$splicing_output}");
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

// second VEP run to annotate RefSeq transcripts
// generate temp file for vep output
$vep_output_refseq = $parser->tempFile("_vep_refseq.vcf");

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


// add phyloP annotation:
$vcf_output_bigwig = $parser->tempFile("_bigwig.vcf");
$parser->exec(get_path("ngs-bits")."/VcfAnnotateFromBigWig", "-name PHYLOP -desc \"".annotation_file_path("/dbs/phyloP/hg38.phyloP100way.bw")." (ngs-bits/VcfAnnotateFromBigWig - mode max)\" -mode max -in {$vep_output_refseq} -out {$vcf_output_bigwig} -bw ".annotation_file_path("/dbs/phyloP/hg38.phyloP100way.bw")." -threads {$threads}", true);


// custom annotation by VcfAnnotateFromVcf

// create config file
$config_file_path = $parser->tempFile(".config");
$config_file = fopen2($config_file_path, 'w');


// add gnomAD annotation
fwrite($config_file, annotation_file_path("/dbs/gnomAD/gnomAD_genome_v3.1.1_GRCh38.vcf.gz")."\tgnomADg\tAF,Hom,Hemi\t\ttrue\n");
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

//perform splicing predictions
if (!$no_splice)
{
	//add missing contig infos
	$start_time = microtime(true);
	add_missing_contigs($vcf_annotate_output);
	$parser->log("Execution time of adding missing contigs: ".time_readable(microtime(true) - $start_time));
	
	//perform splicing annotations
	$start_time = microtime(true);
	annotate_spliceai_score($vcf_annotate_output);
	$parser->log("Execution time of SpliceAI annotation: ".time_readable(microtime(true) - $start_time));
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