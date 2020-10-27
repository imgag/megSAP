<?php 
/** 
	@page vc_freebayes
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_freebayes", "Variant calling with freebayes.");
$parser->addInfileArray("bam",  "Input files in BAM format. Space separated. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output file in VCF.GZ format.", false);
//optional
$parser->addInfile("target",  "Enrichment targets BED file.", true);
$parser->addInt("target_extend",  "Call variants up to n bases outside the target region (they are flagged as 'off-target' in the filter column).", true, 0);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.15);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 1);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 15);
$parser->addInt("min_ao", "Minimum alternative base observation count.", true, 3);
$parser->addFlag("no_ploidy", "Use freebayes parameter -K, i.e. output all alleles which pass input filters, regardles of genotyping outcome or model.");
$parser->addFlag("no_bias", "Use freebayes parameter -V, i.e. ignore strand bias and read end distance bias.");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//create basic variant calls
$args = array();
if(isset($target))
{
	if ($target_extend>0)
	{
		$target_extended = $parser->tempFile("_extended.bed");
		$parser->exec(get_path("ngs-bits")."BedExtend"," -in $target -n $target_extend -out $target_extended -fai {$genome}.fai", true);
	}
	else
	{
		$target_extended = $target;
	}
	
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->exec(get_path("ngs-bits")."BedMerge"," -in $target_extended -out $target_merged", true);
	
	$args[] = "-t $target_merged";
}
if ($no_ploidy)
{
	$args[] = "--pooled-continuous";
}
if ($no_bias)
{
	$args[] = "--binomial-obs-priors-off";
}
$args[] = "--min-alternate-fraction $min_af";
$args[] = "--min-mapping-quality $min_mq";
$args[] = "--min-base-quality $min_bq";
$args[] = "--min-alternate-count $min_ao";
$args[] = "-f $genome";
$args[] = "-b ".implode(" ", $bam);

// run freebayes
$pipeline = array();
if (isset($target) && $threads > 1) 
{	
	$freebayes_start = microtime(true);
	
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
			
			//prepare freebayes arguments
			$args_chr = $args;
			$args_chr[0] = "-t {$tmp_dir}/{$chr}.bed"; //target is always the first parameter
			$args_chr[] = "-v {$tmp_dir}/{$chr}.vcf";
			$stdout_file = "{$tmp_dir}/{$chr}.stdout";
			$args_chr[] = "> $stdout_file";
			$stderr_file = "{$tmp_dir}/{$chr}.stderr";
			$args_chr[] = "2> $stderr_file";
			$parameters = implode(" ", $args_chr);
			
			//start in background
			$command = get_path("freebayes");			
			$output = array();
			$exitcode_file = "{$tmp_dir}/{$chr}.exitcode";
			exec('('.$command.' '.$parameters.' & PID=$!; echo $PID) && (wait $PID; echo $? > '.$exitcode_file.')', $output);
			$pid = trim(implode("", $output));
			
			//log start
			$add_info = array();
			$add_info[] = "version    = ".$parser->extractVersion($command);
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
				$stderr = trim(file_get_contents($stderr_file));
				if ($stderr!="")
				{
					$add_info[] = "STDERR:";
					foreach(explode("\n", $stderr) as $line)
					{
						$add_info[] = nl_trim($line);
					}
				}
				$exit_code = file_get_contents($exitcode_file);
				if($exit_code!=0)
				{
					$add_info[] = "EXIT CODE: ".$exit_code;
				}
				
				//log output
				$parser->log("Finshed processing chromosome {$chr} in ".time_readable(microtime(true)-$start_time), $add_info);
				
				//abort if failed
				if ($exit_code!=0) trigger_error("Processing of chromosome $chr with freebayes failed: ".$stderr, E_USER_ERROR);
				
				unset($running[$chr]);
			}
		}
	}
	$parser->log("Freebayes execution with $threads threads took ".time_readable(microtime(true)-$freebayes_start));
	
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

	$pipeline[] = array("cat", $vcf_combined);
} 
else 
{
	$pipeline[] = array(get_path("freebayes"), implode(" ", $args));
}

//filter variants according to variant quality>5 , alternate observations>=3
$pipeline[] = array(get_path("vcflib")."vcffilter", "-f \"QUAL > 5\"");

//split complex variants to primitives
//this step has to be performed before vcfbreakmulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = array(get_path("vcflib")."vcfallelicprimitives", "-kg");

//split multi-allelic variants
$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");

//normalize all variants and align INDELs to the left
$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize", "-ref $genome");

//sort variants by genomic position
$pipeline[] = array(get_path("ngs-bits")."VcfStreamSort", "");

//fix error in VCF file and strip unneeded information
$pipeline[] = array("php ".repository_basedir()."/src/NGS/vcf_fix.php", "", false);

//zip
$pipeline[] = array("bgzip", "-c > $out", false);

//(2) execute pipeline
$parser->execPipeline($pipeline, "freebayes post processing");

//(3) mark off-target variants
if ($target_extend>0)
{
	$tmp = $parser->tempFile(".vcf");
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $out -mark off-target -reg $target -out $tmp", true);
	$parser->exec("bgzip", "-c $tmp > $out", false);
}

//(4) index output file
$parser->exec("tabix", "-p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return

?>
