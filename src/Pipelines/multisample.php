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
$steps_all = array("vc", "an", "cn");
$parser->addString("steps", "Comma-separated list of steps to perform:\nvc=variant calling, an=annotation, cn=copy-number analysis.", true, implode(",", $steps_all));
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
extract($parser->parse($argv));

//init
$repository_basedir = repository_basedir();
$data_folder = get_path("data_folder");
$hgmd_file = "{$data_folder}/dbs/HGMD/hgmd_cnvs.bed";
$omim_file = "{$data_folder}/dbs/OMIM/omim.bed";

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
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

//(1) variant calling of all samples together (with very conservative parameters)
$vcf_all = $out_folder."all.vcf.gz";
$vcf_all_mito = $out_folder."all_mito.vcf.gz";
$mito = enable_special_mito_vc($sys);
if (in_array("vc", $steps))
{
	$args = array();
	$args[] = "-bam ".implode(" ", $bams);
	$args[] = "-out $vcf_all";
	$args[] = "-target ".$sys['target_file'];
	$args[] = "-min_mq 20";
	$args[] = "-min_af 0.1";
	$args[] = "-target_extend 50";
	$args[] = "-build ".$sys['build'];
	$args[] = "-threads $threads";
	$parser->execTool("NGS/vc_freebayes.php", implode(" ", $args), true);	

	//variant calling for mito with special parameters
	if ($mito)
	{
		$target_mito = $parser->tempFile("_mito.bed");
		file_put_contents($target_mito, "chrMT\t0\t16569");
		
		$args = array();
		$args[] = "-bam ".implode(" ", $bams);
		$args[] = "-out $vcf_all_mito";
		$args[] = "-no_ploidy";
		$args[] = "-min_af 0.01";
		$args[] = "-target $target_mito";
		$args[] = "-build ".$sys['build'];
		$parser->execTool("NGS/vc_freebayes.php", implode(" ", $args), true);
	}
}

//(2) annotation
if (in_array("an", $steps))
{
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
	$parser->exec(get_path("ngs-bits")."VcfStreamSort","-in $vcf -out $vcf_sorted",true);
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
	//check that sample CNV files are present
	$cn_types = array();
	$skip_cn = false;
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
			trigger_error("Skipping CN analysis because sample CNV file is missing: $filename", E_USER_WARNING);
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
	
	//annotation column headers
	$anno_headers = ["overlap af_genomes_imgag", "cn_pathogenic", "genes", "dosage_sensitive_disease_genes", "clinvar_cnvs"];
	
	if (file_exists($hgmd_file))
	{
		$anno_headers[] = "hgmd_cnvs";
	}
	if (file_exists($omim_file))
	{
		$anno_headers[] = "omim";
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
			
			//skip controls;
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
				$match = false;
				foreach($cnv_data[$i] as $cnv2)
				{
					//skip low-confidence CNVs (avoids that noise in unaffected samples removes true-positives in affected samples)
					if ($cnv2["loglike"]<40) continue;
					
					if($cnv["chr"]==$cnv2["chr"] && $cnv["cn"] == $cnv2["cn"] && range_overlap($cnv["start"], $cnv["end"], $cnv2["start"], $cnv2["end"]))
					{
						$match = true;
						break;
					}
				}
				if(!$match)
				{
					$tmp[] = $cnv;
				}
			}
			$regions = $tmp;
		}
		
		$output[] = "#chr\tstart\tend\tsample\tsize\tCN_change\tloglikelihood\tqvalue\tpotential_AF\t".implode("\t", $anno_headers);
		foreach($regions as $reg)
		{
			
			$line = array(
				$reg["chr"],
				$reg["start"],
				$reg["end"],
				$prefix,
				intval($reg["end"])-intval($reg["start"]),
				$reg["cn"],
				$reg["loglike"],
				$reg["qvalue"],
				number_format(max(explode(",", $reg["pot_af"])), 3)
			);
			
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
		$output[] = "#chr	start	end	sample	size	region_count	region_copy_numbers	region_zscores	region_cnv_af	region_coordinates	".implode("\t", $anno_headers);	
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
	
	//(3.2) annotate merged file
	if(!$skip_cn)
	{
		$tmp = $parser->tempFile(".bed");
		file_put_contents($tmp, implode("\n", $output));

		//copy-number polymorphisms
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$tmp} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$tmp}", true);
		
		//knowns pathogenic CNVs
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$tmp} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -out {$tmp}", true);
		
		//gene names
		$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in {$tmp} -extend 20 -out {$tmp}", true);
		
		//dosage sensitive disease genes
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$tmp} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes.bed -no_duplicates -out {$tmp}", true);
		
		//pathogenic ClinVar CNVs
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$tmp} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs.bed -no_duplicates -out {$tmp}", true);
		
		//HGMD CNVs
		if (file_exists($hgmd_file)) //optional because of license
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$tmp} -in2 {$hgmd_file} -no_duplicates -out {$tmp}", true);
		}
		
		//OMIM
		if (file_exists($omim_file)) //optional because of license
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$tmp} -in2 $omim_file -no_duplicates -out {$tmp}", true);
		}
		
		//write output
		if($cn_type == "cnvhunter")	$cnv_multi = "{$out_folder}{$prefix}_cnvs.tsv";
		else if($cn_type == "clincnv") $cnv_multi = "{$out_folder}{$prefix}_cnvs_clincnv.tsv";
		else trigger_error("Invalid CNV list type '{$cn_type}'!", E_USER_ERROR);
		$parser->moveFile($tmp, $cnv_multi);

		//annotate additional gene info
		$parser->exec(get_path("ngs-bits")."CnvGeneAnnotation", "-in {$cnv_multi} -out {$cnv_multi}", true);
	}
}

?>
