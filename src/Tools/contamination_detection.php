<?php
/** 
	@page contamination_detection
	
	Test data: DX180347_01 is contaminated with DX180330_01
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("contamination_detection", "Detects which sample caused the contamination of a second sample.");
$parser->addString("ps", "Processed sample identifier of contaminated sample.", false);
$parser->addEnum("mode", "Select mode to determine variants which are used for comparison.", true, array("AF_outliers_from_GSvar", "AF_outlier_calling"), "AF_outliers_from_GSvar");
$parser->addString("processing_system", "Short name of pocessing system to compare with. If unset, same system as of ps is used.", true);
$parser->addFlag("skip_compare", "Skips comparison to other samples. Use together with 'gsvar_output' option to determine contamination variants.");

// parameters only for mode "AF_outlier_calling"
$parser->addInt("threads", "[mode: AF_outlier_calling] The maximum number of threads used.", true, 1);
$parser->addInfile("system",  "[mode: AF_outlier_calling] Processing system INI file (automatically determined from NGSD if 'ps' is a valid processed sample name).", true);
$parser->addInt("min_var_qual", "[mode: AF_outlier_calling] Minimum quality value for variants.", true,  30);
$parser->addInt("min_var_depth", "[mode: AF_outlier_calling] Minimum depth for variants.", true, 50);
$parser->addInt("min_var_mapq", "[mode: AF_outlier_calling] Minimum mapping quality for variants.", true, 55);
$parser->addFloat("max_af", "[mode: AF_outlier_calling] Maximum gnomAD/1000g allele frequency.", true, 0.05);
$parser->addFloat("max_ngsd", "[mode: AF_outlier_calling] Maximum occurences in NGSD.", true, 100);
$parser->addInfile("ann_vcf", "[mode: AF_outlier_calling] Annotated VCF with variants which should be used for comparison (by a previous run with 'gsvar_output' option).", true);
$parser->addOutfile("gsvar_output", "[mode: AF_outlier_calling] Output GSvar file containing the (filtered) detected variants used for comparison.", true);

extract($parser->parse($argv));

if (is_null($gsvar_output) || $gsvar_output=="")
{
	$gsvar_output = $parser->tempFile(".GSvar");
}

function vcf_lines($filename)
{
	$lines = 0;
	
	$h = fopen($filename, 'r');
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || $line[0]=='#') continue;
		
		++$lines;
	}
	fclose($h);
	
	return $lines;
}

//init
$db = DB::getInstance("NGSD");
$ngsbits = get_path("ngs-bits");
$info = get_processed_sample_info($db, $ps);

if ($mode == "AF_outliers_from_GSvar")
{
	print "##Using mode: AF_outliers_from_GSvar\n";
	//default settings for WES/panel
	$min_dp = 50;
	$r1_start = 0.05;
	$r1_end = 0.30;
	$r2_start = 0.70;
	$r2_end = 0.95;
	if ($info['sys_type']=="WGS")
	{
		print "##WGS sample > adapting parameters\n";
		$min_dp = 30;
		$r1_start = 0.10;
		$r1_end = 0.20;
		$r2_start = 0.80;
		$r2_end = 0.90;
	}
	print "##Using parameters: min_dp={$min_dp} range1={$r1_start}-{$r1_end} range2={$r2_start}-{$r2_end}\n";

	//get variants in unexpected range
	$c_vars = 0;
	$vars = array();
	$gsvar = substr($info['ps_bam'], 0, -4).".GSvar";
	print "##GSvar file: {$gsvar}\n";
	$h = fopen($gsvar, "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || $line[0]=="#") continue;
		
		list($chr, $start, , $ref, $obs, , , $qual) = explode("\t", $line);
		$qual = explode(";", $qual);
		
		//only autosomes
		if (!is_numeric(substr($chr,3))) continue;
		
		//only SNPs
		if (strlen($ref)!=1 || strlen($obs)!=1 || $ref=="-" || $obs=="-") continue;
		
		//only high-depth variants
		$dp = null;
		foreach($qual as $entry)
		{
			if (starts_with($entry, "DP="))
			{
				$dp = substr($entry, 3);
			}
		}
		if ($dp<$min_dp)
		{
			continue;
		}
		
		++$c_vars;
		
		//only strange AF
		$af = null;
		foreach($qual as $entry)
		{
			if (starts_with($entry, "AF="))
			{
				$af = substr($entry, 3);
			}
		}
		if (($af>$r1_start && $af<$r1_end) || ($af>$r2_start && $af<$r2_end))
		{
			$vars["$chr:$start $ref>$obs"] = true;		
		}
	}
	fclose($h);
	print "##Selected ".count($vars)." variants (".number_format(100.0*count($vars)/$c_vars,2)."% of {$c_vars} autosomal SNPs)\n";
}
elseif ($mode == "AF_outlier_calling")
{
	//determine processing system
	$sys = load_system($system, $ps);
	$is_wes = $sys['type']=="WES";
	$is_wgs = $sys['type']=="WGS";
	$is_wgs_shallow = $sys['type']=="WGS (shallow)";
	$has_roi = $sys['target_file']!="";

	//get required files
	$ps_bam = $info['ps_bam'];
	$ps_gsvar = substr($info['ps_bam'], 0, -4).".GSvar";
	$ps_vcf = substr($info['ps_bam'], 0, -4)."_var.vcf.gz";

	// perform additional VC and filtering if no precalculated file is provided
	if (is_null($ann_vcf) || $ann_vcf=="")
	{
		// setup temp files
		$vc_vcf = $parser->tempFile("_vc_noFilter.vcf.gz");
		$decomp_vc_vcf = $parser->tempFile("_vc_decomp.vcf");
		$decomp_ps_vcf = $parser->tempFile("_decomp.vcf");
		$unique_vcf = $parser->tempFile("_unique.vcf");

		//additional VC
		$args = array();
		$args[] = "-bam $ps_bam";
		$args[] = "-out $vc_vcf";
		$args[] = "-threads $threads";
		$args[] = "-min_af 0.01";
		$args[] = "-min_mq 50";
		$args[] = "-min_bq 25";
		$args[] = "-no_ploidy";
		if ($has_roi)
		{
			// limit VC to target region
			$args[] = "-target ".$sys['target_file'];
			$args[] = "-target_extend 50";
		}
		$parser->execTool("NGS/vc_freebayes.php", implode(" ", $args));

		// remove all variants which are also found in original VC

		// decompress original vcf file
		$parser->exec("zcat", "$ps_vcf > $decomp_ps_vcf");

		// decompress new vcf file
		$parser->exec("zcat", "$vc_vcf > $decomp_vc_vcf");

		// parse original VCF and create a array with all variants
		$vars = array();
		$fh = fopen($decomp_ps_vcf, "r");
		$n_var = 0;
		while(!feof($fh))
		{
			$line = trim(fgets($fh));
			// skip comments and header
			if ($line=="" || $line[0]=="#") continue;

			// extract variant info
			list($chr, $pos, , $ref, $obs) = explode("\t", $line);

			//store info:
			$vars["$chr:$pos $ref>$obs"] = true;
			$n_var++;
		}
		fclose($fh);
		print "##Germline VCF file contains {$n_var} variants.\n";

		// parse new vcf and remove all variants which are already found in first run
		$input_fh = fopen($decomp_vc_vcf, "r");
		$output_fh = fopen($unique_vcf, "w");
		$n_var = 0;
		$n_new_var = 0;
		while(!feof($input_fh))
		{
			$line = trim(fgets($input_fh));
			if ($line=="") continue;
			
			// write comments and header
			if ($line[0]=="#") 
			{
				if (starts_with($line, "##"))
				{
					fwrite($output_fh, $line."\n");
				}
				else
				{
					fwrite($output_fh, "##SAMPLE=<ID={$ps},DiseaseStatus=affected>\n");
					fwrite($output_fh, $line."\n");
				}
				continue;
			}
			// extract variant info
			list($chr, $pos, , $ref, $obs) = explode("\t", $line);
			$n_var++;

			// skip all variants which are already found in original VC
			if (isset($vars["$chr:$pos $ref>$obs"])) continue;

			// write unique variants
			fwrite($output_fh, $line."\n");
			$n_new_var++;
		}
		fclose($input_fh);
		fclose($output_fh);
		print "##Low-AF variant calling found $n_var variants.\n";

		// store annotate file in output directory for additional parameter tuning
		$ann_vcf = strtr($gsvar_output, [".GSvar"=>"_ann.vcf"]);
		
		//annotate VCF
		$args = [];
		$args[] = "-in ".$unique_vcf;
		$args[] = "-out ".$ann_vcf;
		$args[] = "-threads ".$threads;
		$args[] = "-ps_name ".$ps;
		$parser->execTool("NGS/an_vep.php", implode(" ", $args));
	}
	else
	{
		$n_new_var = vcf_lines($ann_vcf);
	}
	print "##{$n_new_var} variant found that are not in the germline VCF.\n";

	// filter VCF for variant quality
	$filtered_vcf = $parser->tempFile("_filtered.vcf");
	$parser->exec("{$ngsbits}VcfFilter", "-qual $min_var_qual -sample 'DP >= $min_var_depth' -info 'MQM >= $min_var_mapq' -in $ann_vcf -out $filtered_vcf");
	
	// convert to GSvar
	$gsvar = $parser->tempFile(".GSvar");
	$parser->execTool("NGS/vcf2gsvar.php", "-in $filtered_vcf -out $gsvar");
	print "##".vcf_lines($gsvar)." variant remain after quality filter (QUAL={$min_var_qual}, DP={$min_var_depth}, MQM={$min_var_mapq})\n";
	
	// filter GSvar file for AF and NGSD count
	$filters = [
		"Allele frequency	max_af=".(100.0*$max_af),
		"Allele frequency (sub-populations)	max_af=".(100.0*$max_af),
		"Count NGSD	max_count={$max_ngsd}	ignore_genotype=false",
		"Genotype affected	genotypes=het"
		];
	$filter_file = $parser->tempFile("_filters.txt");
	file_put_contents($filter_file, implode("\n", $filters));
	
	$parser->exec("{$ngsbits}VariantFilterAnnotations", "-filters {$filter_file} -in {$gsvar} -out {$gsvar_output}");
	print "##".vcf_lines($gsvar_output)." variant remain after frequency filter (AF={$max_af}, NGSD={$max_ngsd})\n";
		
	// index filtered variant file
	$vars = array();
	$gsvar_filtered = file($gsvar_output);
	foreach($gsvar_filtered as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		
		// extract variant info
		list($chr, $pos, , $ref, $obs) = explode("\t", $line);

		//store info:
		$vars["$chr:$pos $ref>$obs"] = true;
	}	
}
else
{
	trigger_error("Invalid mode '".$mode."'!", E_USER_ERROR);
}

if ($skip_compare) return;

//determine samples of the same processing system
$tmp = $parser->tempFile(".tsv");
$sys = $info['sys_name_short'];
if(isset($processing_system))
{
	$sys = $processing_system;
}

$pipeline = [
	["{$ngsbits}NGSDExportSamples", "-system {$sys} -no_bad_runs -no_bad_samples -no_tumor -no_ffpe -add_path"],
	["{$ngsbits}TsvFilter", "-filter 'system_name_short is {$sys}'"], //this is necessary because NGSDExportSamples performs fuzzy match
	["{$ngsbits}TsvSlice", "-cols name,path -out {$tmp}"],
];
$parser->execPipeline($pipeline, "NGSD sample extraction", true);
$samples = file($tmp);
print "##Processing ".count($samples)." samples of the processing system {$sys}...\n";

//process samples
print "#sample\toverlap_perc\thom_perc\n";
foreach($samples as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	list($compared_ps, $path) = explode("\t", $line);
	$gsvar = "{$path}/{$compared_ps}.GSvar";
	
	if (!file_exists($gsvar))
	{
		print "##Skipping sample {$compared_ps}: GSvar file missing {$gsvar}\n";
		continue;
	}
	
	//try to re-find variants in the sample
	$c_refound = 0;
	$c_hom = 0;
	$h = fopen($gsvar, "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || $line[0]=="#") continue;
		
		list($chr, $start, , $ref, $obs, $geno) = explode("\t", $line);
		
		if (isset($vars["$chr:$start $ref>$obs"]))
		{
			++$c_refound;
			
			if ($geno=="hom")
			{
				++$c_hom;
			}
		}
	}
	fclose($h);
		
	//print output
	$perc_refound = number_format(100.0*$c_refound/count($vars), 2);
	$hom_perc = $c_refound==0 ? "n/a": number_format(100.0*$c_hom/$c_refound, 2);
	print "{$compared_ps}\t{$perc_refound}\t{$hom_perc}\n";
}

?>
