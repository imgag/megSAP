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
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addString("rna_sample", "Processed sample name of the RNA sample which should be used for annotation.", true, "");
extract($parser->parse($argv));

//input file names
$vcffile = $out_folder."/".$out_name."_var.vcf.gz";
if(!empty($vcf))
{
	$vcffile = $vcf;
}

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
$annfile = $parser->tempFile(".vcf");
$annfile_zipped = $out_folder."/".$out_name."_var_annotated.vcf.gz";	
$varfile = $out_folder."/".$out_name.".GSvar";
$statfile = $out_folder."/".$out_name."_stats_vc.qcML";
$phasing_track = $out_folder."/".$out_name."_phasing_track.bed"; //(only for longreads)


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
$cosmic_cmc = get_path("data_folder") . "/dbs/COSMIC/cmc_export_v98.vcf.gz";
if(file_exists($cosmic_cmc) && $somatic)
{
	$temp_annfile = temp_file(".vcf","cosmic_cmc_an_");
	$parser->exec(get_path("ngs-bits") . "VcfAnnotateFromVcf", "-in {$annfile} -source {$cosmic_cmc} -info_keys COSMIC_CMC -out {$temp_annfile} -threads {$threads}");
	$parser->moveFile($temp_annfile, $annfile);
}

if( $somatic && file_exists(get_path("data_folder")."/dbs/cancerhotspots/cancerhotspots_snv.tsv") )
{
	$temp_annfile = temp_file(".vcf","cosmic_cmc_an_");
	$parser->execTool("NGS/an_somatic_cancerhotspots.php", "-in $annfile -out $temp_annfile");
	$parser->moveFile($temp_annfile, $annfile);
}

//annotate lrGS with short-read WGS sample
if (($sys['type']=="lrGS") && !$multi && !$somatic && db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD");
	//check if out_name is valid ps_name
	if (get_processed_sample_id($db, $out_name, false) != -1)
	{
		list($stdout, $stderr, $exit_code) = $parser->exec(get_path("ngs-bits")."NGSDSameSample", "-ps {$out_name} -system_type WGS");
		//parse stdout
		$same_samples = array();
		foreach ($stdout as $line) 
		{
			if (starts_with($line, "#")) continue;
			if (trim($line) == "") continue;
			$same_samples[] = trim(explode("\t", $line)[0]);
		}
		if (count($same_samples) < 1)
		{
			trigger_error("No related short-read WGS sample found for {$out_name}. Skipping overlap annotation.", E_USER_NOTICE);
		} 
		else
		{
			if (count($same_samples) > 1) trigger_error("Multiple related short-read WGS sample found for {$out_name}. Using first one (".$same_samples[0].").", E_USER_NOTICE);
			else trigger_error("Using short-read WGS sample ".$same_samples[0]." for annotation.", E_USER_NOTICE);
			// get VCF
			list($stdout, $stderr, $exit_code) = $parser->exec(get_path("ngs-bits")."SamplePath", "-ps ".$same_samples[0]." -type VCF", true);
			$sr_vcf = trim($stdout[0]);
			if (!file_exists($sr_vcf)) trigger_error("Short-read VCF '{$sr_vcf}' not found! Skipping annotation.", E_USER_WARNING);
			else
			{
				// annotate variants
				$temp_config = $parser->tempFile("config.tsv");
				file_put_contents($temp_config, "{$sr_vcf}\t\t\t\t\t1\tIN_SHORTREAD_SAMPLE\n");
				$temp_annfile = $parser->tempFile("shortread_ann.vcf");
				$parser->exec(get_path("ngs-bits")."VcfAnnotateFromVcf", "-in {$annfile}  -out {$temp_annfile} -config_file {$temp_config}", true);
				$parser->moveFile($temp_annfile, $annfile);
			}
		}
	}
	else
	{
		trigger_error("Invalid processed sample name provided, cannot get corresponding sr-WGS!", E_USER_WARNING);
	}
	
}


//zip annotated VCF file
$parser->exec("bgzip", "-c $annfile > $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-f -p vcf $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

//convert to GSvar file
if (!$somatic) //germline only
{
	//calculate variant statistics (after annotation because it needs the ID and ANN fields)
	if (!$multi) $parser->exec(get_path("ngs-bits")."VariantQC", "-in $annfile -out $statfile".(($sys['type']=="lrGS")?" -long_read -phasing_bed {$phasing_track}": ""), true);
	
	$args = [];
	$args[] = "-in ".$annfile;
	$args[] = "-out ".$varfile;
	$args[] = "-updown";
	if ($multi) $args[] = "-genotype_mode multi";
	if ($sys['type']=="WGS") $args[] = "-wgs";
	if ($sys['type']=="lrGS") $args[] = "-wgs -longread";
	$parser->execTool("NGS/vcf2gsvar.php", implode(" ", $args));
	
	//annotate variant allele frequency and depth from related RNA sample, if available
	if(db_is_enabled("NGSD"))
	{
		$db = DB::getInstance("NGSD");
		$psample = "";
		if($rna_sample == "")
		{
			$psamples_for_annotation = get_related_processed_samples($db, $out_name, "same sample", "RNA");
			// if multiple RNA samples found: -> choose last one
			if (!empty($psamples_for_annotation)) $psample = end($psamples_for_annotation);
		}
		else
		{
			$psample = $rna_sample;
		}

		if($psample != "")
		{
			trigger_error("Using RNA sample {$psample} for annotation.", E_USER_NOTICE);
			
			$psample_info = get_processed_sample_info($db, $psample);

			//BAM file for variant depth and frequency
			$bam_rna = $psample_info["ps_bam"];
			if (file_exists($bam_rna))
			{
				$parser->exec(get_path("ngs-bits")."VariantAnnotateASE", "-in {$varfile} -out {$varfile} -bam {$bam_rna}", true);

				//check sample similarity
				$min_corr = 0.85;
				$dna_sample_info = get_processed_sample_info($db, $out_name);
				$dna_bam = $dna_sample_info["ps_bam"];

				if (file_exists($dna_bam))
				{
					$output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", "-in {$bam_rna} {$dna_bam} -mode bam -build ".ngsbits_build($sys['build']), true);
					$correlation = explode("\t", $output[0][1])[3];
					if ($correlation=="nan")
					{
						trigger_error("The genotype correlation of DNA and RNA ({$psample}) could not be determined!", E_USER_WARNING);
					}
					else if ($correlation < $min_corr)
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
				$gsvar = Matrix::fromTSV($varfile);
				annotate_gsvar_by_gene($gsvar, $splicing, "symbol", "aberrant_frac", "aberrant_splicing", "Fraction of aberrant splicing reads in gene.");
				$gsvar->toTSV($varfile);
			}
			else
			{
				trigger_error("Splicing file does not exist for RNA sample '{$psample}'!", E_USER_WARNING);
			}

			//add RNA annotation
			$parser->exec(get_path("ngs-bits")."NGSDAnnotateGeneExpression", "-in {$varfile} -out {$varfile} -rna_ps {$psample}", true);
		}
		
	}
	
	//check output TSV file (not for somatic)
	if(!$no_fc)
	{
		$parser->execTool("NGS/check_tsv.php", "-in $varfile -build ".$sys['build']);
	}
}	

?>
