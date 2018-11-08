<?php 
/** 
	@page vc_freebayes
	@todo implement parallelization based on splitting the target region by chromosome
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
$parser->addInt("processes", "How many processes should be used at once to call freebayes", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.15);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 1);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 10);
$parser->addFlag("no_ploidy", "Use freebayes parameter -K, i.e. output all alleles which pass input filters, regardles of genotyping outcome or model.");
$parser->addFlag("no_bias", "Use freebayes parameter -V, i.e. ignore strand bias and read end distance bias.");
extract($parser->parse($argv));

//(1) set up variant calling pipeline
$genome = get_path("local_data")."/{$build}.fa";
$pipeline = array();

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

// run freebayes
if (isset($target) && $processes > 1) 
{	
	// Split BED file by chromosomes into seperate files
	// e.g chr1.bed, chr2.bed, chrY.bed
	$roi = array();
	$bedfile = fopen($target, "r") or die("Cannot read target file: ".$target);
	while (($line = fgets($bedfile)) !== false)
	{
		if (strpos($line, "track") || strpos($line, "browser") || substr_count($line, "\t") > 3) continue;
		$chrom = trim(substr($line, 0, strpos($line, "\t")));
		if (!isset($roi[$chrom])) {
			$roi[$chrom] = array();
		}
		$roi[$chrom][] = $line;
	}
	fclose($bedfile);
	// The list of chromosomes to process for
	$chromosomes = array_keys($roi);

	$tmp_dir = $parser->tempFolder();
	foreach ($chromosomes as $chrom) {
		$stat = file_put_contents($tmp_dir."/".$chrom.".bed", join("", $roi[$chrom]));
		if (!$stat) die("Having trouble saving chromosome ('".$chrom."') to: ".$tmp_dir."/".$chrom);
	}

	/**
	 * Run's the freebayes script with nohup
	 * @function
	 * @param {object} parser
	 * @param {object} bam
	 * @param {object} genome
	 * @param {string} path - the target path to use
	 * @param {string} output - the output path that nohup should redirect stdout to
	 * @return {number} - returns the PID
	 */
	function run_freebayes_nohup($parser, $bam, $genome, $args, $path, $output) {
		$args[0] = "-t ".$path;
		// now we do something like
		// nohup freebayes params &> output &
		// have a look at https://stackoverflow.com/a/4549515/3135319 for further info
		$result = $parser->execParallel("nohup", get_path("freebayes")." -b ".implode(" ",$bam)." -f $genome ".implode(" ", $args)." -v ".$output." </dev/null &", true);
		return $result[3]+1; // returns the PID
	}

	// Then runs the pipeline for every chromosome. Add's n chromosomes to the pool according to the process parameter at the same time.
	$pids = array();
	for ($i = 0; $i < $processes; $i++)
	{
		$chrom = array_shift($chromosomes);
		$pid = run_freebayes_nohup($parser, $bam, $genome, $args, $tmp_dir."/".$chrom.".bed", $tmp_dir."/".$chrom.".vcf");
		array_push($pids, $pid);
	}

	$running = true;
	while ($running) 
	{
		// for all processes check if they are alive
		$running_pids = $parser->exec("pgrep", "-P ".getmypid(), true)[0]; // see https://stackoverflow.com/a/17743940/3135319

		// if all chromosomes have been processed exit the while
		if (!count($chromosomes) && !count($running_pids)) {
			$running = false;
			continue;
		}

		// if less running processes than process limit start a new process
		for ($i = (count($pids) - count($running_pids)); $i < $processes; $i++) 
		{
			if (!count($chromosomes)) continue;
			$chrom = array_shift($chromosomes);
			$pid = run_freebayes_nohup($parser, $bam, $genome, $args, $tmp_dir."/".$chrom.".bed", $tmp_dir."/".$chrom.".vcf");
			array_push($pids, $pid);
		}
	}
	
	// After that merge the resulting VCF files
	$chromosomes = array_keys($roi);
	for ($i = 0; $i < count($chromosomes); $i++) // append all chromsome.vcf files to a combined.vcf
	{
		if ($i != 0) // except for the first chromsome delete all header lines 
		{
			$parser->exec("sed", "-i '/#/d' ".$tmp_dir."/".$chromosomes[$i].".vcf", true);
		}
		$parser->exec("cat", "".$tmp_dir."/".$chromosomes[$i].".vcf >> ".$tmp_dir."/combined.vcf", true);
	}

	unset($roi); // explicitely clean up ROI's because they can be rather large

	// And put a cat on the pipeline script
	$pipeline[] = array("cat", $tmp_dir."/combined.vcf");
} 
else 
{
	$pipeline[] = array(get_path("freebayes"), "-b ".implode(" ",$bam)." -f $genome ".implode(" ", $args));
}

//filter variants according to variant quality>5 , alternate observations>=3
$pipeline[] = array(get_path("vcflib")."vcffilter", "-f \"QUAL > 5 & AO > 2\"");

//split complex variants to primitives
//this step has to be performed before vcfbreakmulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
$pipeline[] = array(get_path("vcflib")."vcfallelicprimitives", "-kg");

//split multi-allelic variants
$pipeline[] = array(get_path("vcflib")."vcfbreakmulti", "");

//normalize all variants and align INDELs to the left
$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize","-ref $genome");

//sort variants by genomic position
$pipeline[] = array(get_path("ngs-bits")."VcfStreamSort","");

//fix error in VCF file and strip unneeded information
$pipeline[] = array("php ".repository_basedir()."/src/NGS/vcf_fix.php", "", false);

//zip
$pipeline[] = array("bgzip", "-c > $out", false);

//(2) execute pipeline
$parser->execPipeline($pipeline, "variant calling");

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
