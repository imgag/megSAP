<?php 
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_varscan2", "Tumor-only variant calling with VarScan2.");
$parser->addInfile("bam",  "Input tumor BAM format.", false);
$parser->addOutfile("out", "Output file in VCF.GZ format (contains both somatic and germline variants).", false);
//optional
$parser->addInfile("target",  "Enrichment target BED file (used to mark off-target calls only).", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFloat("min_af", "Minimum allele frequency cutoff used for variant calling.", true, 0.01);
$parser->addInt("min_dp", "Minimum depth cutoff for variant calling.", true, 20);
$parser->addInt("min_bq", "Minimum base quality cutoff used for variant calling.", true, 20);
$parser->addInt("min_mq", "Minimum mapping quality cutoff used for variant calling.", true, 15);
$parser->addFloat("pval_snv", "p-value threshold for SNVs.", true, 0.05);
$parser->addFloat("pval_indel", "p-value threshold for INDELs.", true, 0.05);
$parser->addString("name", "Sample name to be used in output.", true, "SAMPLE");
$parser->addInfile("debug_region", "Debug option to limit analysis to regions in the given BED file.", true);
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//create mpileup
$pileup_file = $parser->tempFile("_pileup.gz");
$args = [
	"--min-MQ {$min_mq}",
	"--min-BQ {$min_bq}",
	"-f {$genome}",
	"-B", //disable BAQ (increases sensitivitiy significantly
	];
if (isset($debug_region)) $args[] = "-l $debug_region";
$pipeline = [
	["", $parser->execApptainer("samtools", "samtools mpileup", implode(" ", $args)." {$bam}", [$genome, $bam, $debug_region], [], true)],
	["gzip", "--fast > {$pileup_file}"]
];
$parser->execPipeline($pipeline, "mpileup");

//Call SNVs
$tmp_snp_file = $parser->tempFile("varscan2_snps.vcf");
$pipeline = [
	["zcat", $pileup_file],
	["", $parser->execApptainer("varscan2", "java -jar /opt/VarScan.jar", "mpileup2snp --min-var-freq {$min_af} --min-reads2 3 --min_coverage {$min_dp} --output-vcf --p-value {$pval_snv}", [], [], true)],
	["sed", "'/^##source=.*/a ##reference={$genome}'"], //add header line with reference genome after "##source" header line
	["sed", "'s/Sample1/{$name}/'"], //replace sample1 by actual sample $name
	["", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref {$genome} -out {$tmp_snp_file}", [$genome], [], true)]
];
$parser->execPipeline($pipeline, "VarScan2 SNV calling");

//Call INDELs
$tmp_indel_file = $parser->tempFile("varscan2_indels.vcf");
$pipeline = [
	["zcat", $pileup_file],
	["", $parser->execApptainer("varscan2", "java -jar /opt/VarScan.jar", "mpileup2indel --min-var-freq {$min_af} --min-reads2 3 --min_coverage {$min_dp} --output-vcf --p-value {$pval_indel}", [], [], true)],
	["", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref {$genome} -out {$tmp_indel_file}", [$genome], [], true)],
];
$parser->execPipeline($pipeline, "VarScan2 INDEL calling");

//merge SNV and INDEL output
$tmp_merged = $parser->tempFile("varscan2_all_variants.vcf");
$parser->execApptainer("ngs-bits", "VcfAdd", "-in {$tmp_snp_file} {$tmp_indel_file} -out {$tmp_merged}", [] , [dirname($out)]);


//flag variants with AF<5% and calculate QUAL score
$tmp_flagged = $parser->tempFile("varscan2_all_variants.vcf");
$ho = fopen2($tmp_flagged, "w");
$h = fopen2($tmp_merged, "r");
while(!feof($h))
{
	$line = nl_trim(fgets($h));
	if ($line=="") continue;
	
	//headers lines
	if($line[0]=="#")
	{
		if ($line[1]!="#")
		{
			fputs($ho, "##FILTER=<ID=freq-tum,Description=\"Allele frequency less than 5%\">\n");
		}
		fputs($ho, $line."\n");
		continue;
	}
	
	//variant lines
	$parts = explode("\t", $line);
	$format = explode(":", $parts[8]);
	$sample = explode(":", $parts[9]);
	
	//flag low AF flags
	$i_freq = array_search("FREQ", $format);
	if ($i_freq!==false)
	{
		$freq_perc = $sample[$i_freq];
		$freq_perc = substr($freq_perc, 0, -1); //remove % char
		if ($freq_perc<5)
		{
			$filter = $parts[6];
			$parts[6] = ($filter=="PASS" || $filter==".") ? "freq-tum" : $filter.";freq-tum";
		}
	}
	
	//calculate QUAL
	$i_p_value = array_search("PVAL", $format);
	if($i_p_value!==false)
	{
		$snp_q = -10 * log10( $sample[$i_p_value] ); //calculate phred score from p-value
		$snp_q = min(floor($snp_q), 255);
		$parts[5] = $snp_q;
	}
	
	fputs($ho, implode("\t", $parts)."\n");
}
fclose($h);
fclose($ho);

//add version number to source comment:
$vcf = Matrix::fromTSV($tmp_flagged);
$comments = $vcf->getComments();

for ($i=0; $i < count($comments); $i++)
{
	if (contains($comments[$i], "#source="))
	{
		$comments[$i] = $comments[$i]." ".get_path("container_varscan2")."\n";
		break;
	}
}
$comments[] = "#fileDate=".date("Ymd")."\n";

$vcf->setComments($comments);
$vcf->toTSV($tmp_flagged);


//sort output
$parser->execApptainer("ngs-bits", "VcfSort", "-in {$tmp_flagged} -compression_level 9 -out {$out}", [] , [dirname($out)]);

//flag off-targets
if($target!="")
{
	$parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in {$out} -reg {$target} -mark off-target -compression_level 9 -out {$out}", [$out, $target], [dirname($out)]);
}


//index
$parser->execApptainer("htslib", "tabix", "-p vcf $out", [], [dirname($out)]);

?>
