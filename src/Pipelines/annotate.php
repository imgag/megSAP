<?php

/**
	@page annotate
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("annotate", "Annotate variants called on genome build GRCh38.");
$parser->addString("out_name", "Processed sample ID (e.g. 'GS120001_01').", false);
$parser->addString("out_folder", "Output folder.", false);
//optional
$parser->addInfile("system", "Processing system INI file (automatically determined from NGSD if 'out_name' is a valid processed sample name).", true);
$parser->addString("vcf", "Path to (bgzipped) VCF file (if different from {output_folder}/{out_name}_var.vcf.gz).", true, "");
$parser->addString("mosaic_vcf", "Path to (bgzipped) VCF file (if different from {output_folder}/{out_name}_mosaic.vcf.gz).", true, "");
$parser->addFlag("no_fc", "No format check (vcf/tsv).");
$parser->addFlag("multi", "Enable multi-sample mode.");
$parser->addFlag("somatic", "Enable somatic mode (no variant QC and no GSvar file).");
$parser->addFlag("updown", "Don't discard up- or downstream annotations (5000 bases around genes).");
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
extract($parser->parse($argv));


function annotate($parser, $vcffile, $out_name, $out_folder, $system, $no_fc, $multi, $somatic, $updown, $threads, $is_mosaic=false)
{
	//unzip input VCF if necessary
	if (ends_with($vcffile, ".vcf"))
	{
		$vcf_unzipped = $vcffile;
	}
	else
	{
		$vcf_unzipped = $parser->tempFile("_unzipped.vcf");
		$parser->exec("bgzip", "-cd $vcffile > $vcf_unzipped", false); //no output logging, because Toolbase::extractVersion() does not return
	}
	
	//output file names
	$annfile ="";
	$annfile_zipped = "";	
	$varfile = "";
	$stafile = "";
	
	if ($is_mosaic)
	{
		$annfile = $parser->tempFile("mosaic.vcf");
		$annfile_zipped = $out_folder."/".$out_name."_mosaic_annotated.vcf.gz";
		$varfile = $out_folder."/".$out_name."_mosaic.GSvar";
		$stafile = $out_folder."/".$out_name."_mosaic_stats_vc.qcML";
	}
	else
	{
		$annfile = $parser->tempFile(".vcf");
		$annfile_zipped = $out_folder."/".$out_name."_var_annotated.vcf.gz";	
		$varfile = $out_folder."/".$out_name.".GSvar";
		$stafile = $out_folder."/".$out_name."_stats_vc.qcML";
	}

	
	//get system
	$sys = load_system($system, $out_name);
	if ($sys['build']!="GRCh38")
	{
		trigger_error("Unknown genome build ".$sys['build']." cannot be annotated!", E_USER_ERROR);
	}

	//annotate VCF
	$args = [];
	$args[] = "-in ".$vcf_unzipped;
	$args[] = "-out ".$annfile;
	$args[] = "-build ".$sys['build'];
	$args[] = "-threads ".$threads;
	$args[] = "-ps_name ".$out_name;
	if ($somatic) $args[] = "-somatic";
	$parser->execTool("NGS/an_vep.php", implode(" ", $args));

	//annotate COSMIC
	$cosmic_cmc = get_path("data_folder") . "/dbs/COSMIC/cmc_export.vcf.gz";
	if(file_exists($cosmic_cmc) && $somatic)
	{
		$temp_annfile = temp_file(".vcf","cosmic_cmc_an_");
		$parser->exec(get_path("ngs-bits") . "VcfAnnotateFromVcf", "-in {$annfile} -annotation_file {$cosmic_cmc} -info_ids COSMIC_CMC -out {$temp_annfile} -threads {$threads}");
		$parser->moveFile($temp_annfile, $annfile);
	}

	if( $somatic && file_exists(get_path("data_folder")."/dbs/cancerhotspots/cancerhotspots_snv.tsv") )
	{
		$temp_annfile = temp_file(".vcf","cosmic_cmc_an_");
		$parser->execTool("NGS/an_somatic_cancerhotspots.php", "-in $annfile -out $temp_annfile");
		$parser->moveFile($temp_annfile, $annfile);
	}


	//zip annotated VCF file
	$parser->exec("bgzip", "-c $annfile > $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
	$parser->exec("tabix", "-p vcf $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

	//convert to GSvar file
	if (!$somatic) //germline only
	{
		//calculate variant statistics (after annotation because it needs the ID and ANN fields)
		$parser->exec(get_path("ngs-bits")."VariantQC", "-in $annfile -out $stafile", true);
		
		$args = array("-in $annfile", "-out $varfile");
		if ($multi) $args[] = "-genotype_mode multi";
		if ($updown) $args[] = "-updown";
		if ($sys['type']=="WGS") $args[] = "-wgs";
		$parser->execTool("NGS/vcf2gsvar.php", implode(" ", $args));
		
		if(db_is_enabled("NGSD"))
		{
			//annotate variant allele frequency and depth from related RNA sample, if available
			$db = DB::getInstance("NGSD");
			$psamples_for_annotation = get_related_processed_samples($db, $out_name, "same sample", "RNA");
			if (!empty($psamples_for_annotation))
			{
				//last entry
				$psample = end($psamples_for_annotation);
				$psample_info = get_processed_sample_info($db, $psample);

				//BAM file for variant depth and frequency
				$bam_rna = $psample_info["ps_bam"];
				if (file_exists($bam_rna))
				{
					$parser->exec(get_path("ngs-bits")."VariantAnnotateASE", "-in {$varfile} -out {$varfile} -bam {$bam_rna}", true);

					//check sample similarity
					$min_corr = 0.9;
					$dna_sample_info = get_processed_sample_info($db, $out_name);
					$dna_bam = $dna_sample_info["ps_bam"];

					if (file_exists($dna_bam))
					{
						$output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", "-in {$bam_rna} {$dna_bam} -mode bam -build ".ngsbits_build($sys['build']), true);
						$correlation = explode("\t", $output[0][1])[3];
						if ($correlation < $min_corr)
						{
							trigger_error("The genotype correlation of DNA and RNA ({$psample}) is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
						}
						else
						{
							trigger_error("The genotype correlation of DNA and RNA ({$psample}) is {$correlation}.", E_USER_NOTICE);
						}
					}
					else
					{
						trigger_error("BAM file does not exist for DNA sample '{$out_name}'!", E_USER_WARNING);
					}
				}
				else
				{
					trigger_error("BAM file does not exist for RNA sample '{$psample}'!", E_USER_WARNING);
				}

				//load per-gene splicing information
				$splicing = $psample_info['ps_folder']."{$psample}_splicing_gene.tsv";
				if (file_exists($splicing))
				{
					annotate_gsvar_by_gene($varfile, $varfile, $splicing, "symbol", "aberrant_frac", "aberrant_splicing", "Fraction of aberrant splicing reads in gene.");
				}
				else
				{
					trigger_error("Splicing file does not exist for RNA sample '{$psample}'!", E_USER_WARNING);
				}

				//TODO replace with ngs-bits tool:
				//add RNA annotation
				$parser->exec(get_path("ngs-bits")."NGSDAnnotateGeneExpression", "-in {$varfile} -out {$varfile}_ann.GSvar -rna_ps {$psample}", true);

				//load TPM values form counts file, as associative array gene name => TPM value
				$expr = $psample_info['ps_folder']."{$psample}_expr.tsv";
				if (file_exists($expr))
				{
					annotate_gsvar_by_gene($varfile, $varfile, $expr, "gene_name", "tpm", "tpm", "Gene expression strength in transcripts-per-million.");
					annotate_gsvar_by_gene($varfile, $varfile, $expr, "gene_name", "log2fc", "expr_log2fc", "Relative gene expression as log2 FC (log2 tpm).");
					annotate_gsvar_by_gene($varfile, $varfile, $expr, "gene_name", "zscore", "expr_zscore", "Relative gene expression as z-score (log2 tpm)");
				}
				else
				{
					trigger_error("Count file does not exist for RNA sample '{$psample}'!", E_USER_WARNING);
				}

			}
		}
		
		//check output TSV file (not for somatic)
		if(!$no_fc)
		{
			$parser->execTool("NGS/check_tsv.php", "-in $varfile -build ".$sys['build']);
		}
	}	
}

//------------- normal Variants --------------\\

//input file names
$vcffile = $out_folder."/".$out_name."_var.vcf.gz";
if(!empty($vcf))
{
	$vcffile = $vcf;
}

//Annotate the variant file:
annotate($parser, $vcffile, $out_name, $out_folder, $system, $no_fc, $multi, $somatic, $updown, $threads);


//---------------- Mosaic ------------------\\

//annotate mosaic variants file:
$mosaic_file = $out_folder."/".$out_name."_mosaic.vcf.gz";
if (!empty($mosaic_vcf))
{
	$mosaic_file = $mosaic_vcf;
}

if (!file_exists($mosaic_file))
{
	if (! empty($mosaic_vcf))
	{
		trigger_error("Given mosaic vcf file doesn't exist.", E_USER_ERROR);
	}

	return;
}

annotate($parser, $mosaic_file, $out_name, $out_folder, $system, $no_fc, $multi, $somatic, $updown, $threads, true);

?>
