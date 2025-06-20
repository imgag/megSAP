<?php 
/** 
	@page check_for_missing_chromosomes
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("check_for_missing_chromosomes", "Check for missing chromosomes (chr1-chr22 and chrX only).");
$parser->addInfile("in",  "Input VCF or VCF.GZ file.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFloat("max_missing_perc", "If set, the percentage of bases without variants for each chromosome is checked (used for small variants mainly).", true, 0);
$parser->addFlag("debug", "Add debug output");
extract($parser->parse($argv));

//init
$chrs = chr_list(false);
$is_gvcf = ends_with(strtolower($in), ".gvcf") || ends_with(strtolower($in), ".gvcf.gz");

//determine variant coordiante range for each chromosome
$ranges = [];
$h = gzopen2($in, "r");
while(!gzeof($h))
{
	$line = trim(gzgets($h));
	if ($line=="" || $line[0]=="#") continue;
	
	$parts = explode("\t", $line, 6);
	$chr = $parts[0];
	$pos = $parts[1];
	if($chr=="chrMT" || $chr=="chrY") continue;
	if ($is_gvcf && $parts[4]=="<NON_REF>") continue;
	
	if (!isset($ranges[$chr])) $ranges[$chr] = [$pos, $pos];
	if ($pos<$ranges[$chr][0]) $ranges[$chr][0] = $pos;
	if ($ranges[$chr][1]<$pos) $ranges[$chr][1] = $pos;
}
gzclose($h);

if ($debug)
{
	foreach($ranges as $chr => list($from, $to))
	{
		print "DEBUG - coverered: $chr:$from-$to\n";
	}
}

//check for missing chromosomes
$missing = [];
foreach ($chrs as $chr) 
{
	if (!isset($ranges[$chr])) $missing[] = $chr;
}
if (count($missing)>0) trigger_error("Missing chromosomes ".implode(", ", $missing)." in {$in}", E_USER_ERROR);

//check for chromosome parts without variants
if ($max_missing_perc>0)
{
	//get chromosome sizes
	$sizes = [];
	list($stdout) = exec2("cut -f1,2 ".genome_fasta($build).".fai");
	foreach($stdout as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		if (!in_array($chr, $chrs)) continue;
		
		list($chr, $size) = explode("\t", $line);
		$sizes[$chr] = $size;
	}

	//check covered chromosome range
	$missing = [];
	foreach($ranges as $chr => list($from, $to))
	{
		//calculae missing bases
		$tmp_covered = temp_file(".bed");
		file_put_contents($tmp_covered, "{$chr}\t{$from}\t{$to}");
		$tmp_missing = temp_file(".bed");
		list($stdout) = exec2("echo '{$chr}\t1\t".$sizes[$chr]."' | BedSubtract -in2 {$tmp_covered} | BedSubtract -in2 ".repository_basedir()."/data/misc/nobase_regions.bed | BedInfo | grep Bases");
		$bases_missing = explode(":", $stdout[0])[1];
		
		$bases_missing_perc = number_format(100.0 * $bases_missing / $sizes[$chr], 1);
		if ($debug) print "DEBUG - missing bases: $chr:$from-$to missing=$bases_missing_perc%\n";
		if ($bases_missing_perc>$max_missing_perc) $missing[] = "{$chr}={$bases_missing_perc}%";
	}

	if (count($missing)>0) trigger_error("Chromosome percentage without variants: ".implode(", ", $missing)." in {$in}", E_USER_ERROR);

}

?>
