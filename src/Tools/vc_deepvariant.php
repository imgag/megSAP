<?php 
/** 
	@page vc_deepvariant
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_deepvariant", "Variant calling with DeepVariant.");
$parser->addInfileArray("bam",  "Input files in BAM format. Space separated. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output file in VCF.gz format.", false);
$parser->addString("model_type", "Type of model to use for variant calling. Choose from <WGS|WES|PACBIO|ONT_R104|HYBRID_PACBIO_ILLUMINA|MASSEQ>.", false);
//optional
$parser->addInfile("target",  "Enrichment targets BED file.", true);
$parser->addInt("target_extend",  "Call variants up to n bases outside the target region (they are flagged as 'off-target' in the filter column).", true, 0);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("gvcf", "Enable output of gVCF files and define output filepath for gVCF files.", true, "");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.15);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 20);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 10);
$parser->addFlag("raw_output", "return the raw output of deepvariant with no post-processing.");
$parser->addFlag("allow_empty_examples", "allows DeepVariant to call variants even if no examples were created with make_examples.");
$parser->addFlag("gpu", "Use GPU supportet DeepVariant container");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//create basic variant calls
$args = array();
$in_files = array();
if(isset($target))
{	
	//extend by 'n' bases
	$target_extended = $parser->tempFile("_roi_extended.bed");
	if ($target_extend>0)
	{
		$parser->execApptainer("ngs-bits", "BedExtend", "-in $target -n $target_extend -out $target_extended -fai {$genome}.fai", [$target, $genome]);
	}
	else
	{
		$parser->copyFile($target, $target_extended);
	}
	
	//add special target regions (regions with known pathogenic variants that are often captured by exome/panel, but not inside the target region)
	if ($build=="GRCh38" && $target_extend>0) //only if extended (otherwise it is also added for chrMT calling, etc.)
	{
		$parser->execApptainer("ngs-bits", "BedAdd", "-in $target_extended ".repository_basedir()."data/misc/special_regions.bed -out $target_extended ", [repository_basedir()."data/misc/special_regions.bed"]);
	}
	
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->execApptainer("ngs-bits", "BedMerge", "-in $target_extended -out $target_merged");
	
	$args[] = "--regions $target_merged";
}
if ($allow_empty_examples)
{
	$args[] = "--call_variants_extra_args=allow_empty_examples=true";
}

$args[] = "--model_type=$model_type";
$args[] = "--make_examples_extra_args=min_mapping_quality=$min_mq,min_base_quality=$min_bq,vsc_min_fraction_indels=$min_af,vsc_min_fraction_snps=$min_af";
$args[] = "--ref=$genome";
$args[] = "--reads=".implode(" ", $bam);
$args[] = "--num_shards=".$threads;

if (!empty($gvcf)) $args[] = "--output_gvcf=$gvcf";

$in_files[] = $genome;
$in_files = array_merge($in_files, $bam);

// run deepvariant
$pipeline = array();
$container = ($gpu) ? "deepvariant-gpu" : "deepvariant";
if (isset($target) && $threads > 20) //TODO set back to > 1 if prefered over num_shards
{	
	$deepvar_start = microtime(true);
	
	// split BED file by chromosomes into seperate files  e.g chr1.bed, chr2.bed, chrY.bed
	$roi_by_chr = array();
	$chr_order_original = array();
	$chr_order_by_size = array();
	$file = file($target_merged);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || starts_with($line, "track") || starts_with($line, "browser") || substr_count($line, "\t") < 2) continue;
		list($chr, $start, $end) = explode("\t", $line);
		if (!isset($roi_by_chr[$chr]))
		{
			$roi_by_chr[$chr] = array();
			$chr_order_original[] = $chr;
			$chr_order_by_size[$chr] = 0;
		}
		$roi_by_chr[$chr][] = $line."\n";
		$chr_order_by_size[$chr] += $end - $start;
	}
	
	// store each chromosome in separate file
	$tmp_dir = $parser->tempFolder("vc_freebayes_pid".getmypid()."_"); //we identify sub-processes using this folder name > include the PID
	foreach ($roi_by_chr as $chr => $lines)
	{
		file_put_contents("{$tmp_dir}/{$chr}.bed", $lines);
	}
	unset($roi_by_chr);
	
	//sort chromosomes by size
	arsort($chr_order_by_size);
	$chr_order_by_size = array_keys($chr_order_by_size);
	
	// run variant calling for every chromosome separately
	$running = array();
	$chrs = $chr_order_by_size;
	while (true) 
	{
		//wait a second
		sleep(1);
		
		//all chromosomes have been processed > exit
		if (count($chrs)==0 && count($running)==0) break;
		
		//start new sub-processed while there are threads unused
		while(count($running)<$threads && count($chrs)>0)
		{
			$chr = array_shift($chrs);
			
			//prepare deepvariant arguments
			$args_chr = $args;
			$args_chr[0] = "--regions={$tmp_dir}/{$chr}.bed"; //target is always the first parameter
			$args_chr[] = "--output_vcf={$tmp_dir}/{$chr}.vcf";
			$stdout_file = "{$tmp_dir}/{$chr}.stdout";
			$args_chr[] = "> $stdout_file";
			$stderr_file = "{$tmp_dir}/{$chr}.stderr";
			$args_chr[] = "2> $stderr_file";
			$parameters = implode(" ", $args_chr);
			
			//start in background
			$command = $parser->execApptainer($container, "run_deepvariant", $parameters, $in_files, [], true);		
			$output = array();
			$exitcode_file = "{$tmp_dir}/{$chr}.exitcode";
			exec('('.$command.' & PID=$!; echo $PID) && (wait $PID; echo $? > '.$exitcode_file.')', $output);
			$pid = trim(implode("", $output));
			
			//log start
			$add_info = array();
			$add_info[] = "version    = ".get_path("container_deepvariant");
			$add_info[] = "parameters = {$parameters}";
			$add_info[] = "pid = {$pid}";
			$parser->log("Executing command in background '{$command}'", $add_info);
			
			$running[$chr] = array($pid, $stdout_file, $stderr_file, $exitcode_file, microtime(true));
		}
		
		//check which started processes are still runnning
		$chrs_running = array_keys($running);
		foreach($chrs_running as $chr)
		{
			list($pid, $stdout_file, $stderr_file, $exitcode_file, $start_time) = $running[$chr];
			if(!file_exists("/proc/{$pid}"))
			{
				//prepare additional infos for logging
				$add_info = array();
				$stdout = trim(file_get_contents($stdout_file));
				if ($stdout!="")
				{
					$add_info[] = "STDOUT:";
					foreach(explode("\n", $stdout) as $line)
					{
						$add_info[] = nl_trim($line);
					}
				}
				$job_aborted = false;
				$stderr = trim(file_get_contents($stderr_file));
				if ($stderr!="")
				{
					$add_info[] = "STDERR:";
					foreach(explode("\n", $stderr) as $line)
					{
						$add_info[] = nl_trim($line);
						if (contains($line, "terminate called after throwing an instance of 'std::runtime_error'"))
						{
							$job_aborted = true;
						}
					}
				}
				$exit_code = trim(file_get_contents($exitcode_file));
				$add_info[] = "EXIT CODE: ".$exit_code;
				if(!is_numeric($exit_code) || $exit_code!=0)
				{
					$job_aborted = true;
				}
				
				//log output
				$parser->log("Finshed processing chromosome {$chr} in ".time_readable(microtime(true)-$start_time), $add_info);
				
				//abort if failed
				if ($job_aborted) trigger_error("Processing of chromosome $chr with deepvariant failed: ".$stderr, E_USER_ERROR);
				
				unset($running[$chr]);
			}
		}
	}
	$parser->log("Deepvariant execution with $threads threads took ".time_readable(microtime(true)-$deepvar_start));
	
	// combine individual chromosome VCFs to one file
	$combine_start = microtime(true);
	$vcf_combined = "{$tmp_dir}/combined.vcf";
	$ho = fopen2($vcf_combined, "w");
	for ($i = 0; $i < count($chr_order_original); $i++)
	{
		$h = fopen2("{$tmp_dir}/".$chr_order_original[$i].".vcf", "r");
		while(!feof($h))
		{
			$line = fgets($h);
			if (trim($line)=="") continue;
			
			//skip headers (except for first chromosome)
			if ($line[0]=='#')
			{
				if ($i==0)
				{
					fwrite($ho, $line);
				}
				continue;
			}
			
			fwrite($ho, $line);
		}
		fclose($h);
	}
	fclose($ho);
	$parser->log("Combining VCFs took ".time_readable(microtime(true)-$combine_start));
	
	// if only raw output is requested copy combined_vcf to output and stop:
	if ($raw_output)
	{
		exec("cp $vcf_combined $out");
		return;
	}
	
	$pipeline[] = array("cat", $vcf_combined);

	//filter variants according to variant quality>5
	$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfFilter", "-qual 5 -remove_invalid -ref $genome", [$genome], [], true)];
} 
else 
{
	if ($raw_output)
	{
		$parser->execApptainer($container, "run_deepvariant" ,implode(" ", $args)." --output_vcf=$out", $in_files, [dirname($out)]);
		return;
	}

	$vcf_deepvar_out = $parser->tempFile(".vcf");
	$parser->execApptainer($container, "run_deepvariant", implode(" ", $args)." --output_vcf=$vcf_deepvar_out", $in_files, [dirname($out)]);

	//filter variants according to variant quality>5
	$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfFilter", "-in $vcf_deepvar_out -qual 5 -remove_invalid -ref $genome", [$genome], [], true)];
}



//split complex variants to primitives
//this step has to be performed before VcfBreakMulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = ["", $parser->execApptainer("vcflib", "vcfallelicprimitives", "-kg", [], [], true)];

//split multi-allelic variants - -no_errors flag can be removed, when vcfallelicprimitives is replaced
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "-no_errors", [], [], true)];
//normalize all variants and align INDELs to the left
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref $genome", [$genome], [], true)];

//sort variants by genomic position
$tmp_out = $parser->tempFile(".vcf");
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfStreamSort", "-out $tmp_out", [], [], true)];

//(2) execute pipeline
$parser->execPipeline($pipeline, "deepvariant post processing");

//prepend source, date and reference to outfile:
$file_format = "##fileformat=VCFv4.2\n";
$file_date = "##fileDate=".date("Ymd")."\n";
$source_line = "##source=DeepVariant ".get_path("container_deepvariant")."\n";
$reference_line = "##reference=".genome_fasta($build, false)."\n";
file_put_contents($tmp_out, $file_format . $file_date . $source_line . $reference_line . file_get_contents($tmp_out));

//zip
$parser->exec("bgzip", "-c $tmp_out > $out");

//(3) mark off-target variants
if ($target_extend>0)
{
	$tmp = $parser->tempFile(".vcf");
	$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in $out -mark off-target -reg $target -out $tmp", [$out, $target]);
	$parser->exec("bgzip", "-c $tmp > $out", false);
}

//(4) index output file
$parser->exec("tabix", "-f -p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return

?>
