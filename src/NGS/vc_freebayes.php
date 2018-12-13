<?php 
/** 
	@page vc_freebayes
	@todo test if '--use-best-n-alleles 4' speed up processing and if results are still ok.
	@todo test hard filters SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1, see https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=6&ved=2ahUKEwj4k4zdrvndAhUuM-wKHaCXC40QFjAFegQIAxAC&url=https%3A%2F%2Fwiki.uiowa.edu%2Fdownload%2Fattachments%2F145192256%2Ferik%2520garrison%2520-%2520iowa%2520talk%25202.pdf%3Fapi%3Dv2&usg=AOvVaw0G6VgcVVuS42Bk2WBlP1IS
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
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 10);
$parser->addFlag("no_ploidy", "Use freebayes parameter -K, i.e. output all alleles which pass input filters, regardles of genotyping outcome or model.");
$parser->addFlag("no_bias", "Use freebayes parameter -V, i.e. ignore strand bias and read end distance bias.");
extract($parser->parse($argv));

//init
$genome = get_path("local_data")."/{$build}.fa";

//create basic variant calls
$args = array();
if(isset($target))
{
	if ($target_extend>0)
	{
		$target_extended = $parser->tempFile("_extended.bed");
		$parser->exec(get_path("ngs-bits")."BedExtend"," -in $target -n $target_extend -out $target_extended -fai ".get_path("data_folder")."/genomes/".$build.".fa.fai", true);
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
$args[] = "--min-base-quality $min_bq"; //max 10% error propbability
$args[] = "--min-alternate-qsum 90"; //At least 3 good observations
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
			$args_chr = $args;
			$args_chr[0] = "-t {$tmp_dir}/{$chr}.bed"; //target is always the first parameter
			$args_chr[] = "-v {$tmp_dir}/{$chr}.vcf";
			list($stdout, $stderr, $status) = $parser->execBackground(get_path("freebayes"), implode(" ", $args_chr));
			$running[$chr] = array($status['pid'], $stdout, $stderr, microtime(true));
		}
		
		//check which started processes are still runnning
		$tmp_base = basename($tmp_dir);
		list($processes) = exec2("ps ax | grep ".get_path("freebayes"));
		$chrs_running = array_keys($running);
		foreach($chrs_running as $chr)
		{
			$still_running = false;
			foreach($processes as $process)
			{
				if (contains($process, $tmp_base) && contains($process, "/{$chr}.vcf"))
				{
					$still_running = true;
				}
			}
			if(!$still_running)
			{
				//check for error
				$stderr = trim(file_get_contents($running[$chr][2]));
				if ($stderr!="")
				{
					if(stripos($stderr, "error")!==FALSE)
					{
						trigger_error("Processing of chromosome $chr with freebayes failed: ".$stderr, E_USER_ERROR);
					}
					else
					{
						trigger_error("Processing of chromosome $chr returned the following output on STDERR: ".$stderr, E_USER_WARNING);
					}
				}
				
				//print execution time
				$parser->log("Finshed processing chromosome {$chr} in ".time_readable(microtime(true)-$running[$chr][3]));
				
				unset($running[$chr]);
			}
		}
	}
	$parser->log("Freebayes execution with $threads threads took ".time_readable(microtime(true)-$freebayes_start));
	
	// combine individual chromosome VCFs to one file
	$combine_start = microtime(true);
	$vcf_combined = "{$tmp_dir}/combined.vcf";
	$ho = fopen($vcf_combined, "w");
	for ($i = 0; $i < count($chr_order_original); $i++)
	{
		$h = fopen("{$tmp_dir}/".$chr_order_original[$i].".vcf", "r");
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

$freebayes_end = microtime(true);
$parser->log("Freebayes execution took ".time_readable($freebayes_end-$freebayes_start));

//filter variants according to variant quality>5 (~31% error probabilty) , alternate allele observations>2
$pipeline[] = array(get_path("ngs-bits")."VcfFilter", "-qual 5 -info \"AO > 2\"");

//split multi-allelic variants
$pipeline[] = array(get_path("ngs-bits")."VcfBreakMulti", "");

//split complex variants to primitives
$pipeline[] = array(get_path("ngs-bits")."VcfBreakComplex", "");

//normalize all variants and align INDELs to the left
$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize","-ref $genome");

//sort variants by genomic position
$pipeline[] = array(get_path("ngs-bits")."VcfStreamSort","");

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
