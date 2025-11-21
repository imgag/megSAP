##fileformat=VCFv4.2
##INFO=<ID=ALPHAMISSENSE,Number=.,Type=String,Description="AlphaMissense score from AlphaMissense_hg38.tsv.gz">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$h = gzopen($argv[1], 'r');
while (!feof($h))
{
	$line = trim(fgets($h));
	$parts = explode("\t", $line);
	if ($line=="" || $line[0]=="#") continue;
	if (count($parts)==10) //for AlphaMissense_hg38.tsv.gz
	{
		list($chr, $pos, $ref, $alt, $build, $prot, $enst, $aa, $score, $class) = $parts;
	}
	else if (count($parts)==9) //for AlphaMissense_isoforms_hg38.vcf.gz
	{
		list($chr, $pos, $ref, $alt, $build, $enst, $aa, $score, $class) = $parts;
	}
	else trigger_error("Invalid column count in line $line", E_USER_ERROR);
	
	if ($chr=="chrM") $chr="chrMT";
	print implode("\t", [$chr, $pos, ".", $ref, $alt, "30", "PASS", "ALPHAMISSENSE=$score"])."\n";
}
gzclose($h);

?>