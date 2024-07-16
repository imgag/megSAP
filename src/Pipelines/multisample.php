<?php

/**
	@page multisample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function extract_info($format, $data)
{
	if ($data==".")
	{
		return array('NONE', 0, 0);
	}
	
	$data = explode(":", $data);
	$data = array_combine($format, $data);
	$depth = array_sum(explode(",", $data['DP'])); 
	$genotype = vcfgeno2human($data['GT']);
	$ao = array_sum(explode(",", $data['AO']));
	return array($genotype, $depth, $ao);
}


//parse command line arguments
$parser = new ToolBase("multisample", "Multi-sample analysis pipeline.");
$parser->addInfileArray("bams", "Input BAM files.", false);
$parser->addStringArray("status", "List of affected status of the input samples (BAMs) - can be 'affected' or 'control'.", false);
$parser->addString("out_folder", "Output folder name.", false);
//optional
$parser->addString("prefix", "Output file prefix.", true, "multi");
$parser->addInfile("system",  "Processing system INI file used for all samples (automatically determined from NGSD if the basename of the first file in 'bams' is a valid processed sample name).", true);
$steps_all = array("vc", "cn", "sv", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nvc=variant calling, cn=copy-number analysis, sv=structural variant calling, db=database import.", true, implode(",", $steps_all));
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("annotation_only", "Performs only a reannotation of the already created variant calls.");
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
extract($parser->parse($argv));

//init
$repository_basedir = repository_basedir();
$data_folder = get_path("data_folder");
$ngsbits = get_path("ngs-bits");
$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2023_3.bed";
$omim_file = "{$data_folder}/dbs/OMIM/omim.bed";
$vcf_all = "{$out_folder}/all.vcf.gz";
$vcf_all_mito = "{$out_folder}/all_mito.vcf.gz";
$cnv_multi = "{$out_folder}/{$prefix}_cnvs_clincnv.tsv";
$gsvar = "{$out_folder}/{$prefix}.GSvar";
$sv_manta_file = "{$out_folder}/{$prefix}_manta_var_structural.vcf.gz";
$bedpe_out = substr($sv_manta_file,0,-6)."bedpe";

// create logfile in output folder if no filepath is provided:
if ($parser->getLogFile() == "") $parser->setLogFile($out_folder."/multi_".date("YmdHis").".log");

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//log server name
list($server) = exec2("hostname -f");
$user = exec('whoami');
$parser->log("Executed on server: ".implode(" ", $server)." as ".$user);

if ($annotation_only)
{
	// check if required VC files for annotation is available and print warning otherwise
	if (in_array("vc", $steps) && !file_exists($vcf_all))
	{
		trigger_error("VCF for reannotation is missing. Skipping 'vc' step!", E_USER_WARNING);
		if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
	}

	if (in_array("cn", $steps) && !file_exists($cnv_multi))
	{
		$cnv_multi = "{$out_folder}{$prefix}_cnvs_clincnv.tsv";
		if (!file_exists($cnv_multi))
		{
			trigger_error("CN file for reannotation is missing. Skipping 'cn' step!", E_USER_WARNING);
			if (($key = array_search("cn", $steps)) !== false) unset($steps[$key]);
		}
	}

	if (in_array("sv", $steps) && !file_exists($bedpe_out))
	{
		trigger_error("BEDPE file for reannotation is missing. Skipping 'sv' step!", E_USER_WARNING);
		if (($key = array_search("sv", $steps)) !== false) unset($steps[$key]);
	}
}

//check input counts
if (count($bams)!=count($status))
{
	trigger_error("Different number of arguments for 'bams' and 'status' parameters! BAMs: ".count($bams)." Status: ".count($status)."!", E_USER_ERROR);
}

//check input status
foreach($status as $stat)
{
	$valid = array("affected", "control");
	if (!in_array($stat, $valid))
	{
		trigger_error("Invalid status '$stat' given for '$bam'. Valid are '".implode("', '", $valid)."'", E_USER_ERROR);
	}
}
$status = array_combine($bams, $status);

//check input sample names
$names = array();
$tmp = array();
foreach($bams as $bam)
{
	$name = basename2($bam);
	if (isset($tmp[$name]))
	{
		trigger_error("Sample file name '$name' occurs twice in input file list. Each sample must be uniquely indentifyable by name!", E_USER_ERROR);
	}
	$names[$bam] = $name;
	$tmp[$name] = true;
}


//check processing systems are valid and matching
$sys = load_system($system, $names[$bams[0]]);
$sys_matching = true;
foreach($bams as $bam)
{
	$sys_ps = load_system($system, $names[$bam]);
	if ($sys_ps['target_file']=="")
	{
		trigger_error("Cannot perform multi-sample analysis without target region (processing systems of ".$bam." is '".$sys["name_short"]."')!", E_USER_ERROR);
	}
	if ($sys["name_short"]!=$sys_ps["name_short"]) $sys_matching = false;
}

//check genome build of BAMs
foreach($bams as $bam)
{
	check_genome_build($bam, $sys['build']);
}


//check steps
$is_wgs = $sys['type']=="WGS";
$is_wgs_shallow = $sys['type']=="WGS (shallow)";
if ($is_wgs_shallow)
{
	if (in_array("vc", $steps))
	{
		trigger_error("Skipping step 'vc' - Variant calling is not supported for shallow WGS samples!", E_USER_NOTICE);
		if (($key = array_search("vc", $steps)) !== false) unset($steps[$key]);
	}
	if (in_array("an", $steps))
	{
		trigger_error("Skipping step 'an' - Annotation is not supported for shallow WGS samples!", E_USER_NOTICE);
		if (($key = array_search("an", $steps)) !== false) unset($steps[$key]);
	}
}

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync)
{
	$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);
}

//create output folder if missing
$out_folder .= "/";
if (!file_exists($out_folder))
{
	if (!mkdir($out_folder))
	{
		trigger_error("Could not create output folder '{$out_folder}'!", E_USER_ERROR);
	}
	if (!chmod($out_folder, 0777))
	{
		trigger_error("Could not change privileges of output folder '{$out_folder}'!", E_USER_ERROR);
	}
}

//copy BAM files to local tmp for variant calling
if (!$annotation_only && (in_array("vc", $steps) || in_array("sv", $steps)))
{
	$tmp_bam_folder = $parser->tempFolder("local_copy_of_bams_");
	$local_bams = array();
	foreach($bams as $bam)
	{
		//convert to BAM if CRAM
		$bam = convert_to_bam_if_cram($bam, $parser, $sys['build'], $threads, $tmp_bam_folder);
		
		//check BAM/BAI exist
		if (!file_exists($bam)) trigger_error("BAM file '{$bam}' does not exist!", E_USER_ERROR);
		if (!file_exists($bam.".bai")) trigger_error("BAM index file '{$bam}.bai' does not exist!", E_USER_ERROR);
		
		//create local copy if not already local file
		if (starts_with($bam, $tmp_bam_folder))
		{
			$local_bams[] = $bam;
		}
		else
		{
			$local_bam = $tmp_bam_folder."/".basename($bam);
			$parser->copyFile($bam, $local_bam);
			$parser->copyFile($bam.".bai", $local_bam.".bai");
			$local_bams[] = $local_bam;
		}
	}
}

//(1) variant calling of all samples together (with very conservative parameters)
$mito = enable_special_mito_vc($sys);
if (in_array("vc", $steps))
{
	if (!$annotation_only)
	{
		$args = array();
		$args[] = "-bam ".implode(" ", $local_bams);
		$args[] = "-out $vcf_all";
		$args[] = "-target ".$sys['target_file'];
		$args[] = "-min_mq 20";
		$args[] = "-min_af 0.1";
		$args[] = "-target_extend 200";
		$args[] = "-build ".$sys['build'];
		$args[] = "-threads $threads";
		$args[] = "--log ".$parser->getLogFile();
		$parser->execTool("NGS/vc_freebayes.php", implode(" ", $args), true);	

		//variant calling for mito with special parameters
		if ($mito)
		{
			$target_mito = $parser->tempFile("_mito.bed");
			file_put_contents($target_mito, "chrMT\t0\t16569");
			
			$args = array();
			$args[] = "-bam ".implode(" ", $local_bams);
			$args[] = "-out $vcf_all_mito";
			$args[] = "-no_ploidy";
			$args[] = "-no_bias";
			$args[] = "-min_af 0.01";
			$args[] = "-target $target_mito";
			$args[] = "-build ".$sys['build'];
			$args[] = "--log ".$parser->getLogFile();
			$parser->execTool("NGS/vc_freebayes.php", implode(" ", $args), true);
		}
	}
	//(2) annotation
	//Convert VCF to single-sample format
	$indices = array();
	$h1 = gzopen2($vcf_all, "r");
	$vcf = $parser->tempFile("_unzipped.vcf");
	$h2 = fopen2($vcf, "w");
	while(!gzeof($h1))
	{
		$line = trim(gzgets($h1));
		if (strlen($line)==0) continue;
		
		if ($line[0]=="#" && $line[1]=="#") //comments
		{
			fwrite($h2, $line."\n");
		}
		else if ($line[0]=="#") //header
		{
			//add multi-sample comments
			fwrite($h2, "##FORMAT=<ID=MULTI,Number=.,Type=String,Description=\"Multi-sample genotype information (genotype, depth, alternate base observations).\">\n");		
			fwrite($h2, "##ANALYSISTYPE=GERMLINE_MULTISAMPLE\n");
			fwrite($h2, "##PIPELINE=".repository_revision(true)."\n");
			foreach($bams as $bam)
			{
				fwrite($h2, gsvar_sample_header($names[$bam], array("DiseaseStatus"=>$status[$bam])));
			}
			
			//determine indices for each sample	
			$parts = explode("\t", $line);
			foreach($bams as $bam)
			{
				$indices[$bam] = vcf_column_index($names[$bam], $parts); 
			}
			
			//write main header
			fwrite($h2, implode("\t", array_slice($parts, 0, 9))."\t{$prefix}\n");
		}
		else //content
		{
			$parts = explode("\t", $line);
			$format = explode(":", $parts[8]);
			$muti_info = array();
			foreach($bams as $bam)
			{
				$index = $indices[$bam];
				list($gt, $dp, $ao) = extract_info($format, $parts[$index]);
				$muti_info[] = $names[$bam]."=$gt|$dp|$ao";
			}
			
			//update format field
			$parts[8] = "MULTI";
			fwrite($h2, implode("\t", array_slice($parts, 0, 9))."\t".implode(",", $muti_info)."\n");
		}
	}
	if ($mito && file_exists($vcf_all_mito)) //check if file exists for downward compatibility with old analyses that did not include mito
	{
		//clean up
		fclose($h1);
		$indices = array();
		$h1 = gzopen2($vcf_all_mito, "r");
		while(!gzeof($h1))
		{
			$line = trim(gzgets($h1));
			if (strlen($line)==0) continue;
			
			if ($line[0]=="#") //header > determine indices
			{
				if ($line[1]=="#") continue; //comments
				
				$parts = explode("\t", $line);
				foreach($bams as $bam)
				{
					$indices[$bam] = vcf_column_index($names[$bam], $parts); 
				}
			}
			else //content > append
			{
				$parts = explode("\t", $line);
				$format = explode(":", $parts[8]);
				$muti_info = array();
				foreach($bams as $bam)
				{
					$index = $indices[$bam];
					list($gt, $dp, $ao) = extract_info($format, $parts[$index]);
					$muti_info[] = $names[$bam]."=$gt|$dp|$ao";
				}
				
				//update format field
				$parts[8] = "MULTI";
				fwrite($h2, implode("\t", array_slice($parts, 0, 9))."\t".implode(",", $muti_info)."\n");
			}
		}
	}
	fclose($h1);
	fclose($h2);

	//Sort variant list (neccessary for tabix)
	$vcf_sorted = $parser->tempFile("_unsorted.vcf");
	$parser->exec("{$ngsbits}VcfStreamSort","-in $vcf -out $vcf_sorted",true);
	//zip variant list
	$vcf_zipped = "{$out_folder}{$prefix}_var.vcf.gz";
	$parser->exec("bgzip", "-c $vcf_sorted > $vcf_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
	$parser->exec("tabix", "-p vcf $vcf_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

	//basic annotation
	$parser->execTool("Pipelines/annotate.php", "-out_name {$prefix} -out_folder $out_folder -system $system -threads $threads -multi");

	//update sample entry 
	$status_map = array();
	foreach ($status as $bam => $disease_status) 
	{
		$status_map[basename2($bam)] = $disease_status;
	}
	update_gsvar_sample_header($gsvar, $status_map);

}

//(3) copy-number calling
if (in_array("cn", $steps))
{
	//(3.1) CNV calling
	if (!$annotation_only)
	{
		//load target region BED for exome (to determine the number of regions)
		$roi = [];
		if (!$is_wgs && !$is_wgs_shallow && $sys_matching)
		{
			foreach(file($sys['target_file']) as $line)
			{
				$line = trim($line);
				if ($line=="" || $line[0]=="#") continue;
				list($chr, $start, $end) = explode("\t", $line);
				$roi[$chr][] = [$start, $end];
			}
		}
		
		//(3.1a) merge ClinCNV calls of individual samples
		$output = array();
		$output[] = "##ANALYSISTYPE=CLINCNV_GERMLINE_MULTI";
		$cnv_data = array();
		foreach($bams as $bam)
		{
			$ps_name = basename2($bam);
			
			//parse CNV file
			$data = array();
			$cnv_file = dirname($bam)."/{$ps_name}_cnvs_clincnv.tsv";
			if (!file_exists($cnv_file)) trigger_error("CNV file missing: {$cnv_file}", E_USER_ERROR);
			foreach(file($cnv_file) as $line)
			{
				$line = nl_trim($line);
				if($line == "") continue;
				
				//header
				if(starts_with($line, "#"))
				{
					if(starts_with($line,"##"))
					{
						if (starts_with($line, "##ANALYSISTYPE=")) continue;
						
						$output[] = "##{$ps_name} ".substr($line, 2);
					}
					continue;
				}
				
				//content line
				$parts = explode("\t",$line);
				list($chr,$start,$end,$cn,$loglikelihood,$no_of_regions,$length_KB,$potential_af,$genes,$qvalue) = explode("\t", $line);
				
				$data[] = array("chr" => $chr, "start" => $start, "end" =>$end, "cn" => $cn, "loglike" =>$loglikelihood, "pot_af" => $potential_af, "qvalue" => $qvalue);
			}
			$cnv_data[] = $data;
		}
		
		//intersect CNV regions of affected (with same CN state)
		$regions = null;
		$i = -1;
		foreach($status as $bam => $stat)
		{
			++$i;
			
			//skip controls
			if ($stat != "affected") continue;
			
			//first affected => init
			if (is_null($regions))
			{
				$regions = $cnv_data[$i];
				continue;
			}
			
			$tmp = array();
			foreach($regions as $cnv)
			{
				foreach($cnv_data[$i] as $cnv2)
				{
					if($cnv["chr"]==$cnv2["chr"] && $cnv["cn"] == $cnv2["cn"] && range_overlap($cnv["start"]+1, $cnv["end"], $cnv2["start"]+1, $cnv2["end"]))
					{
						$new_entry = $cnv;
						$new_entry["start"] = max($cnv["start"], $cnv2["start"]);
						$new_entry["end"] = min($cnv["end"], $cnv2["end"]);
						$new_entry["loglike"] .= ",".$cnv2["loglike"];
						$new_entry["qvalue"] .= ",".$cnv2["qvalue"];
						$new_entry["pot_af"] = max($cnv["pot_af"], $cnv2["pot_af"]);
						
						$tmp[] = $new_entry;
						break;
					}
				}
			}
			$regions = $tmp;
		}
		
		//subtract CNV regions of controls (with same CN state)
		$i=-1;
		foreach($status as $bam => $stat)
		{
			++$i;
			
			//skip affected
			if($stat == "affected") continue;
			
			//subtract
			$tmp = array();
			foreach($regions as $cnv)
			{
				$match_bases = 0;
				foreach($cnv_data[$i] as $cnv2)
				{
					//skip low-confidence CNVs (avoids that noise in unaffected samples removes true-positives in affected samples)
					if ($cnv2["loglike"]<20) continue;
					
					if($cnv["chr"]==$cnv2["chr"] && $cnv["cn"] == $cnv2["cn"] && range_overlap($cnv["start"], $cnv["end"], $cnv2["start"], $cnv2["end"]))
					{
						$overlap_start = max($cnv["start"], $cnv2["start"]);
						$overlap_end = min($cnv["end"], $cnv2["end"]);
						$match_bases += $overlap_end - $overlap_start;
						break;
					}
				}
				$min_match_bases = 0.5 * ($cnv["end"] - $cnv["start"]);
				if($match_bases < $min_match_bases)
				{
					$tmp[] = $cnv;
				}
			}
			$regions = $tmp;
		}
		
		$output[] = "#".implode("\t", ["chr", "start", "end", "sample", "size", "CN_change", "loglikelihood", "no_of_regions", "qvalue", "potential_AF"]);
		foreach($regions as $reg)
		{
			$chr = $reg["chr"];
			$start = $reg["start"];
			$end = $reg["end"];
			$size = intval($end)-intval($start);
			$line = array(
				$chr,
				$start,
				$end,
				$prefix,
				$size,
				$reg["cn"],
				$reg["loglike"]
				);
			if($is_wgs)
			{
				$bin_size = get_path("cnv_bin_size_wgs");
				$line[] = round($size/$bin_size); 
			}
			else if ($is_wgs_shallow)
			{
				$bin_size = get_path("cnv_bin_size_shallow_wgs");
				$line[] = round($size/$bin_size); 
			}
			else if ($sys_matching) //calculate based on the number of target regions overlapping
			{
				$ol_count = 0;
				if (isset($roi[$chr]))
				{
					foreach($roi[$chr] as list($start2, $end2))
					{
						if (range_overlap($start, $end, $start2+1, $end2-1)) ++$ol_count;
					}
				}
				$line[] = $ol_count; 
			}
			else
			{
				$line[] = "n/a"; 
			}
			
			$line[] = $reg["qvalue"];
			$line[] = number_format(max(explode(",", $reg["pot_af"])), 3);
			
			$output[] = implode("\t",$line);
		}
	
		//write output file
		file_put_contents($cnv_multi, implode("\n", $output));

		//create dummy GSvar file for shallow WGS (needed to be able to open the sample in GSvar)
		if ($is_wgs_shallow)
		{
			$content = array();
			$content[] = "##ANALYSISTYPE=GERMLINE_TRIO";
			foreach($bams as $bam)
			{
				$content[] = trim(gsvar_sample_header($names[$bam], array("DiseaseStatus"=>$status[$bam])));
			}
			$desc_and_filter = array(
				"##DESCRIPTION=filter=Annotations for filtering and ranking variants.",
				"##DESCRIPTION=quality=Quality parameters - Quality parameters - variant quality (QUAL), depth (DP), allele frequency (AF), mean mapping quality of alternate allele (MQM), probability of strand bias for alternate bases as phred score (SAP), probability of allele ballance as phred score (ABP).",
				"##DESCRIPTION=gene=Affected gene list (comma-separated).",
				"##DESCRIPTION=variant_type=Variant type.",
				"##DESCRIPTION=coding_and_splicing=Coding and splicing details (Gene, ENST number, type, impact, exon/intron number, HGVS.c, HGVS.p, Pfam domain).",
				"##DESCRIPTION=regulatory=Regulatory consequence details.",
				"##DESCRIPTION=OMIM=OMIM database annotation.",
				"##DESCRIPTION=ClinVar=ClinVar database annotation.",
				"##DESCRIPTION=HGMD=HGMD database annotation.",
				"##DESCRIPTION=RepeatMasker=RepeatMasker annotation.",
				"##DESCRIPTION=dbSNP=Identifier in dbSNP database.",
				"##DESCRIPTION=1000g=Allele frequency in 1000 genomes project.",
				"##DESCRIPTION=gnomAD=Allele frequency in gnomAD project.",
				"##DESCRIPTION=gnomAD_hom_hemi=Homoyzgous counts and hemizygous counts of gnomAD project (genome data).",
				"##DESCRIPTION=gnomAD_sub=Sub-population allele frequenciens (AFR,AMR,EAS,NFE,SAS) in gnomAD project.",
				"##DESCRIPTION=phyloP=phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6.",
				"##DESCRIPTION=CADD=CADD pathogenicity prediction scores (scaled phred-like). Deleterious threshold > 10-20.",
				"##DESCRIPTION=REVEL=REVEL pathogenicity prediction score. Deleterious threshold > 0.5.",
				"##DESCRIPTION=AlphaMissense=AlphaMissense pathogenicity prediction score. Deleterious threshold > 0.564.",
				"##DESCRIPTION=MaxEntScan=MaxEntScan splicing prediction (reference bases score/alternate bases score).",
				"##DESCRIPTION=COSMIC=COSMIC somatic variant database anntotation.",
				"##DESCRIPTION=NGSD_hom=Homozygous variant count in NGSD.",
				"##DESCRIPTION=NGSD_het=Heterozygous variant count in NGSD.",
				"##DESCRIPTION=NGSD_group=Homozygous / heterozygous variant count in NGSD with the same disease group (Neoplasms).",
				"##DESCRIPTION=classification=Classification from the NGSD.",
				"##DESCRIPTION=classification_comment=Classification comment from the NGSD.",
				"##DESCRIPTION=validation=Validation information from the NGSD. Validation results of other samples are listed in brackets!",
				"##DESCRIPTION=comment=Variant comments from the NGSD.",
				"##DESCRIPTION=gene_info=Gene information from NGSD (inheritance mode, gnomAD o/e scores).",
				"##FILTER=low_conf_region=Low confidence region for small variant calling based on gnomAD AC0/RF filters and IMGAG trio/twin data.",
				"##FILTER=off-target=Variant marked as 'off-target'."
				);
			$content = array_merge($content, $desc_and_filter);
			$content[] = "#chr	start	end	ref	obs	".implode("\t", array_values($names))."	filter	quality	gene	variant_type	coding_and_splicing	regulatory	OMIM	ClinVar	HGMD	RepeatMasker	dbSNP	1000g	gnomAD	gnomAD_hom_hemi	gnomAD_sub	phyloP	CADD	REVEL	MaxEntScan	COSMIC	NGSD_hom	NGSD_het	NGSD_group	classification	classification_comment	validation	comment	gene_info";
			file_put_contents($gsvar, implode("\n", $content));
		}
	}

	
	//(3.2) annotate merged file
	//copy-number polymorphisms
	$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$cnv_multi}", true);
	
	//knowns pathogenic CNVs
	$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -out {$cnv_multi}", true);

	//annotate additional gene info
	$parser->exec("{$ngsbits}CnvGeneAnnotation", "-in {$cnv_multi} -add_simple_gene_names -out {$cnv_multi}", true);
	
	//dosage sensitive disease genes
	$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed -no_duplicates -out {$cnv_multi}", true);
	
	//pathogenic ClinVar CNVs
	$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2024-02.bed -name clinvar_cnvs -url_decode -no_duplicates -out {$cnv_multi}", true);
	
	//HGMD CNVs
	if (file_exists($hgmd_file)) //optional because of license
	{
		$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$cnv_multi}", true);
	}
	
	//OMIM
	if (file_exists($omim_file)) //optional because of license
	{
		$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 $omim_file -no_duplicates -url_decode -out {$cnv_multi}", true);
	}
}

//sv calling
if (in_array("sv", $steps))
{
	// skip SV calling if only annotation should be done	
	if (!$annotation_only)
	{
		//SV calling with manta
		$manta_evidence_dir = "{$out_folder}/manta_evid";
		create_directory($manta_evidence_dir);


		$manta_args = [
			"-bam ".implode(" ", $local_bams),
			"-evid_dir ".$manta_evidence_dir,
			"-out ".$sv_manta_file,
			"-threads ".$threads,
			"-build ".$sys['build'],
			"--log ".$parser->getLogFile()
		];
		$manta_args[] = "-target ".$sys['target_file'];
		if(!$is_wgs) $manta_args[] = "-exome";
		
		$parser->execTool("NGS/vc_manta.php", implode(" ", $manta_args));

		// Rename Manta evidence file
		$evidence_bam_files = glob("$manta_evidence_dir/evidence_*.bam");
		foreach($evidence_bam_files as $old_bam_filename)
		{
			$new_bam_filename = preg_replace(array("/evidence_[0-9]+\./", "/\.bam/"), array("", "_manta_evidence.bam" ), $old_bam_filename);
			$parser->moveFile($old_bam_filename, $new_bam_filename);
			$parser->moveFile($old_bam_filename.".bai", $new_bam_filename.".bai");
		}
		
		//create BEDPE files
		$parser->exec("{$ngsbits}VcfToBedpe", "-in $sv_manta_file -out $bedpe_out", true);

		// correct filetype and add sample headers
		$bedpe_table = Matrix::fromTSV($bedpe_out);
		$bedpe_table->removeComment("#fileformat=BEDPE");
		if ($prefix == "trio")
		{
			$bedpe_table->prependComment("#fileformat=BEDPE_GERMLINE_TRIO");
		}
		else
		{
			$bedpe_table->prependComment("#fileformat=BEDPE_GERMLINE_MULTI");
		}
		$bedpe_table->addComment("#PIPELINE=".repository_revision(true));
		foreach($bams as $bam)
		{
			$bedpe_table->addComment(gsvar_sample_header($names[$bam], array("DiseaseStatus"=>$status[$bam]), "#"));
		}
		$bedpe_table->toTSV($bedpe_out);	
	}

	//add gene info annotation and NGSD counts
	if (db_is_enabled("NGSD"))
	{
		$parser->exec("{$ngsbits}BedpeGeneAnnotation", "-in $bedpe_out -out $bedpe_out -add_simple_gene_names", true);
	}

	//add NGSD counts from flat file
	$ngsd_annotation_folder = get_path("data_folder")."/dbs/NGSD/";
	$ngsd_sv_files = array("sv_deletion.bedpe.gz", "sv_duplication.bedpe.gz", "sv_insertion.bedpe.gz", "sv_inversion.bedpe.gz", "sv_translocation.bedpe.gz");
	$db_file_dates = array();

	// check file existance
	$all_files_available = file_exists($ngsd_annotation_folder."sv_breakpoint_density.igv");
	foreach ($ngsd_sv_files as $filename) 
	{
		if(!(file_exists($ngsd_annotation_folder.$filename) && file_exists($ngsd_annotation_folder.$filename.".tbi")))
		{
			$all_files_available = false;
			break;
		}
	}
	if ($all_files_available)
	{
		// store flat file modification date to detect changes during annotation 
		foreach ($ngsd_sv_files as $filename)
		{
			$db_file_dates[$filename] = filemtime($ngsd_annotation_folder.$filename);
			if ($db_file_dates[$filename] == false)
			{
				trigger_error("Cannot get modification date of '".$ngsd_annotation_folder.$filename."'!",E_USER_ERROR);
			}
		}
		
		//perform annotation
		$parser->exec("{$ngsbits}BedpeAnnotateCounts", "-in $bedpe_out -out $bedpe_out -processing_system ".$sys["name_short"]." -ann_folder {$ngsd_annotation_folder}", true);
		$sys_specific_density_file = $ngsd_annotation_folder."sv_breakpoint_density_".$sys["name_short"].".igv";
		if (file_exists($sys_specific_density_file))
		{
			$parser->exec("{$ngsbits}BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv -density_sys {$sys_specific_density_file}", true);
		}
		else
		{
			$parser->exec("{$ngsbits}BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv", true);
		}

		// check if files changed during annotation
		foreach ($ngsd_sv_files as $filename)
		{
			if ($db_file_dates[$filename] != filemtime($ngsd_annotation_folder.$filename))
			{
				trigger_error("Annotation file '".$ngsd_annotation_folder.$filename."' has changed during annotation!",E_USER_ERROR);
			}
		}

	}
	else
	{
		trigger_error("Cannot annotate NGSD counts! At least one required file in '{$ngsd_annotation_folder}' is missing!", E_USER_WARNING);
	}

	//add optional OMIM annotation

	$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; 
	if(file_exists($omim_file))//OMIM annotation (optional because of license)
	{
		$parser->exec("{$ngsbits}BedpeAnnotateFromBed", "-in $bedpe_out -out $bedpe_out -bed $omim_file -url_decode -replace_underscore -col_name OMIM", true);
	}

	//add CNV overlap annotation
	if (file_exists($cnv_multi))
	{
		$parser->exec("{$ngsbits}BedpeAnnotateCnvOverlap", "-in $bedpe_out -out $bedpe_out -cnv $cnv_multi", true);
	}

	//update sample entry 
	$status_map = array();
	foreach ($status as $bam => $disease_status) 
	{
		$status_map[basename2($bam)] = $disease_status;
	}
	update_gsvar_sample_header($bedpe_out, $status_map);
}

//NGSD import
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	//add secondary analysis (if missing)
	$parser->execTool("NGS/db_import_secondary_analysis.php", "-type 'multi sample' -gsvar {$gsvar}");
}

?>
