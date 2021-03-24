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
extract($parser->parse($argv));

//init
$repository_basedir = repository_basedir();
$data_folder = get_path("data_folder");
$ngsbits = get_path("ngs-bits");
$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2020_4.bed";
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
		$cnv_multi = "{$out_folder}{$prefix}_cnvs.tsv";
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
	$name = basename($bam, ".bam");
	if (isset($tmp[$name]))
	{
		trigger_error("Sample file name '$name' occurs twice in input file list. Each sample must be uniquely indentifyable by name!", E_USER_ERROR);
	}
	$names[$bam] = $name;
	$tmp[$name] = true;
}

//check processing systems are supported
foreach($bams as $bam)
{
	$sys = load_system($system, $names[$bam]);
	if ($sys['target_file']=="")
	{
		trigger_error("Cannot perform multi-sample analysis without target region (processing systems of ".$bam." is '".$sys["name_short"]."')!", E_USER_ERROR);
	}
}


$sys = load_system($system, $names[$bams[0]]);

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

//create output folder if missing
$out_folder .= "/";
if (!file_exists($out_folder))
{
	mkdir($out_folder);
	if (!chmod($out_folder, 0777))
	{
		trigger_error("Could not change privileges of folder '{$out_folder}'!", E_USER_ERROR);
	}
}

// copy BAM files to local tmp 
if (!$annotation_only && (in_array("vc", $steps) || in_array("sv", $steps)))
{
	$local_bams = array();
	$tmp_bam_folder = $parser->tempFolder("bam_");
	foreach($bams as $bam)
	{
		if (!file_exists($bam)) trigger_error("BAM file '$bam' not found in Sample folder! Cannot perform any calling steps!", E_USER_ERROR);
		if (!file_exists($bam.".bai")) trigger_error("BAM index file for BAM '$bam' not found in Sample folder!", E_USER_ERROR);
		$local_bam = $tmp_bam_folder."/".basename($bam);
		$parser->copyFile($bam, $local_bam);
		$parser->copyFile($bam.".bai", $local_bam.".bai");
		$local_bams[] = $local_bam;
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
		$args[] = "-target_extend 50";
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

}

//(3) copy-number calling
if (in_array("cn", $steps))
{
	$skip_cn = false;
	if (!$annotation_only)
	{		
		//check that sample CNV files are present
		$cn_types = array();
		foreach($bams as $bam)
		{
			$base = substr($bam, 0, -4);
			$filename_cnvhunter = $base."_cnvs.tsv";
			$filename_clincnv = $base."_cnvs_clincnv.tsv";
			
			if (file_exists($filename_clincnv))
			{
				$cn_types[] = "clincnv";
			}
			else if (file_exists($filename_cnvhunter))
			{
				$cn_types[] = "cnvhunter";
			}
			else
			{
				trigger_error("Skipping CN analysis because sample CNV file is missing: $filename_clincnv", E_USER_WARNING);
				$skip_cn = true;
				break;
			}
		}
		
		//check that all CNV files are of the same type
		$cn_types = array_unique($cn_types);
		if(count($cn_types) == 1)
		{
			$cn_type = $cn_types[0];
		}
		else
		{
			trigger_error("Input CNVs have different filetypes. Aborting CNV multisample call.",E_USER_WARNING);
			$skip_cn = true;
		}
		
		if(!$skip_cn)
		{
			if($cn_type == "cnvhunter")	$cnv_multi = "{$out_folder}/{$prefix}_cnvs.tsv";
			else if($cn_type == "clincnv") $cnv_multi = "{$out_folder}/{$prefix}_cnvs_clincnv.tsv";
			else trigger_error("Invalid CNV list type '{$cn_type}'!", E_USER_ERROR);
		}
		
		//(3.1a) merge ClinCNV
		if(!$skip_cn && $cn_type == "clincnv")
		{
			//load CNV data
			$output = array();
			$output[] = "##ANALYSISTYPE=CLINCNV_GERMLINE_MULTI";
			$cnv_data = array();
			foreach($bams as $bam)
			{
				$ps_name = basename($bam, ".bam");
				
				//parse CNV file
				$data = array();
				$file = file(substr($bam, 0, -4)."_cnvs_clincnv.tsv");
				foreach($file as $line)
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
			
			$headers = array("chr", "start", "end", "sample", "size", "CN_change", "loglikelihood");
			if($is_wgs || $is_wgs_shallow) $headers[] = "no_of_regions";
			$headers[] = "qvalue";
			$headers[] = "potential_AF";
			$output[] = "#".implode("\t", $headers);
			foreach($regions as $reg)
			{
				$size = intval($reg["end"])-intval($reg["start"]);
				$line = array(
					$reg["chr"],
					$reg["start"],
					$reg["end"],
					$prefix,
					$size,
					$reg["cn"],
					$reg["loglike"]
					);
				if($is_wgs || $is_wgs_shallow)
				{
					$bin_size = get_path("cnv_bin_size_wgs");
					if ($is_wgs_shallow) $bin_size = get_path("cnv_bin_size_shallow_wgs");
					$line[] = round($size/$bin_size); 
				}
				$line[] = $reg["qvalue"];
				$line[] = number_format(max(explode(",", $reg["pot_af"])), 3);
				
				$output[] = implode("\t",$line);
			}
		}
		
		//(3.1b) merge CnvHunter
		if (!$skip_cn && $cn_type == "cnvhunter")
		{
			//load CNV files
			$output = array();
			$output[] = "##ANALYSISTYPE=CNVHUNTER_GERMLINE_MULTI";
			$cnv_data = array();
			foreach($bams as $bam)
			{
				//parse CNV file
				$cnv_num = 0;
				$data = array();
				$file = file(substr($bam, 0, -4)."_cnvs.tsv");
				foreach($file as $line)
				{
					$line = nl_trim($line);
					if ($line=="") continue;
					
					//header line
					if (starts_with($line, "#"))
					{
						if (starts_with($line, "##"))
						{
							if (starts_with($line, "##ANALYSISTYPE=")) continue;
							
							$output[] = $line;
						}
						continue;
					}
					
					//content line
					$parts = explode("\t", $line);
					$regs = explode(",", $parts[9]);
					$cns = explode(",", $parts[6]);
					$zs = explode(",", $parts[7]);
					$afs = explode(",", $parts[8]);
					for($i=0; $i<count($regs); ++$i)
					{
						$data[] = array($regs[$i], $cns[$i], $zs[$i], $afs[$i], $cnv_num);
					}
					++$cnv_num;
				}
				$cnv_data[] = $data;
			}
			
			//intersect CNV regions of affected (with same CN state)
			$regions = null;
			$i=-1;
			foreach($status as $bam => $stat)
			{
				++$i;
				
				//skip controls
				if ($stat!="affected") continue;
				
				//first affected => init
				if (is_null($regions))
				{
					$regions = $cnv_data[$i];
					continue;
				}
				
				//intersect
				$tmp = array();
				foreach($regions as $cnv)
				{
					foreach($cnv_data[$i] as $cnv2)
					{
						if ($cnv[0]==$cnv2[0] && $cnv[1]==$cnv2[1])
						{
							$new_entry = $cnv;
							$new_entry[2] .= ",".$cnv2[2];
							$tmp[] = $new_entry;
							break;
						}
					}
				}
				$regions = $tmp;
			}
			
			//subract CNV regions of controls (with same CN state)
			$i=-1;
			foreach($status as $bam => $stat)
			{
				++$i;
				
				//skip affected
				if ($stat=="affected") continue;
						
				//subtract
				$tmp = array();
				foreach($regions as $cnv)
				{
					$match = false;
					foreach($cnv_data[$i] as $cnv2)
					{
						if ($cnv[0]==$cnv2[0] && $cnv[1]==$cnv2[1])
						{
							$match = true;
							break;
						}
					}
					if (!$match)
					{
						$tmp[] = $cnv;
					}
				}
				$regions = $tmp;
			}
			
			//re-group regions by original number ($cnv_num)
			$regions_by_number = array();
			foreach($regions as $cnv)
			{
				$regions_by_number[$cnv[4]][] = $cnv;
			}
			
			//write output
			$output[] = "#chr	start	end	sample	size	region_count	region_copy_numbers	region_zscores	region_cnv_af	region_coordinates";	
			foreach($regions_by_number as $cnv_num => $regions)
			{
				$chr = null;
				$start = null;
				$end = null;
				$cns = array();
				$zs = array();
				$afs = array();
				$cords = array();
				foreach($regions as $reg)
				{
					list($c, $s, $e) = explode(":", strtr($reg[0], "-", ":"));
					
					//coordinate range
					if (is_null($chr))
					{
						$chr = $c;
						$start = $s;
						$end = $e;
					}
					else
					{
						$start = min($start, $s);
						$end = max($end, $e);
					}
					
					$cns[] = $reg[1];
					$z_scores = explode(",", $reg[2]);
					$zs[] = mean($z_scores);
					$afs[] = $reg[3];
					$cords[] = $reg[0];
					
				}
				$output[] = "{$chr}\t{$start}\t{$end}\t{$prefix}\t".($end-$start)."\t".count($regions)."\t".implode(",", $cns)."\t".implode(",", $zs)."\t".implode(",", $afs)."\t".implode(",", $cords);
			}
		}

		//write output file
		if (!$skip_cn)
		{
			file_put_contents($cnv_multi, implode("\n", $output));
		}

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
				"##DESCRIPTION=Sift=Sift effect prediction and score for each transcript: D=damaging, T=tolerated.",
				"##DESCRIPTION=PolyPhen=PolyPhen (humVar) effect prediction and score for each transcript: D=probably damaging, P=possibly damaging, B=benign.",
				"##DESCRIPTION=fathmm-MKL=fathmm-MKL score (for coding/non-coding regions). Deleterious threshold > 0.5.",
				"##DESCRIPTION=CADD=CADD pathogenicity prediction scores (scaled phred-like). Deleterious threshold > 10-20.",
				"##DESCRIPTION=REVEL=REVEL pathogenicity prediction score. Deleterious threshold > 0.5.",
				"##DESCRIPTION=MaxEntScan=MaxEntScan splicing prediction (reference bases score/alternate bases score).",
				"##DESCRIPTION=dbscSNV=dbscSNV splicing prediction (ADA/RF score).",
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
			$content[] = "#chr	start	end	ref	obs	".implode("\t", array_values($names))."	filter	quality	gene	variant_type	coding_and_splicing	regulatory	OMIM	ClinVar	HGMD	RepeatMasker	dbSNP	1000g	gnomAD	gnomAD_hom_hemi	gnomAD_sub	phyloP	Sift	PolyPhen	fathmm-MKL	CADD	REVEL	MaxEntScan	dbscSNV	COSMIC	NGSD_hom	NGSD_het	NGSD_group	classification	classification_comment	validation	comment	gene_info";
			file_put_contents($gsvar, implode("\n", $content));
		}
	}

	
	//(3.2) annotate merged file
	if(!$skip_cn)
	{
		//copy-number polymorphisms
		$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$cnv_multi}", true);
		
		//knowns pathogenic CNVs
		$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -out {$cnv_multi}", true);

		//annotate additional gene info
		$parser->exec("{$ngsbits}CnvGeneAnnotation", "-in {$cnv_multi} -out {$cnv_multi} -add_simple_gene_names", true);
		
		//dosage sensitive disease genes
		$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes.bed -no_duplicates -out {$cnv_multi}", true);
		
		//pathogenic ClinVar CNVs
		$parser->exec("{$ngsbits}BedAnnotateFromBed", "-in {$cnv_multi} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2021-01.bed -name clinvar_cnvs -url_decode -no_duplicates -out {$cnv_multi}", true);
		
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
		if($sys['target_file']!="") $manta_args[] = "-target ".$sys['target_file'];
		if(!$is_wgs) $manta_args[] = "-exome";
		
		$parser->execTool("NGS/vc_manta.php", implode(" ", $manta_args));

		// Rename Manta evidence file
		$evidence_bam_files = glob("$manta_evidence_dir/evidence_*.bam");
		foreach($evidence_bam_files as $old_bam_filename)
		{
			$new_bam_filename = preg_replace(array("/evidence_[0-9]+\./", "/\.bam/"), array("", "_manta_evidence.bam" ), $old_bam_filename);
			rename($old_bam_filename, $new_bam_filename);
			rename($old_bam_filename.".bai", $new_bam_filename.".bai");
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
		$parser->exec("{$ngsbits}NGSDAnnotateSV", "-in $bedpe_out -out $bedpe_out -ps $name", true);
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
}

//NGSD import
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	//add secondary analysis (if missing)
	$parser->execTool("NGS/db_import_secondary_analysis.php", "-type 'multi sample' -gsvar {$gsvar}");
}

?>
