<?php
/** 
	@page kasp_ngs_match 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//returns the genotype(s) for a sample at a certain position, or 'n/a' if the minimum depth was not reached.
function ngs_geno($bam, $chr, $pos, $ref, $min_depth)
{
	//check if cached
	static $cache = [];
	$key = $chr . ":" . $pos . " " . $ref;
	if(isset($cache[$bam][$key]))
	{
		return $cache[$bam][$key];
	}
	
	//get pileup
	list($output) = exec2(get_path("samtools")." mpileup -aa -r $chr:$pos-$pos $bam");
	list($chr2, $pos2, $ref2, , $bases) = explode("\t", $output[0]);;
	
	//count bases
	$bases = strtoupper($bases);
	$counts = array("A"=>0, "C"=>0, "G"=>0, "T"=>0);
	for($i=0; $i<strlen($bases); ++$i)
	{
		$char = $bases[$i];
		if (isset($counts[$char]))
		{
			++$counts[$char];
		}
	}
	arsort($counts);
	
	//check depth
	$keys = array_keys($counts);
	$c1 = $counts[$keys[0]];
	$c2 = $counts[$keys[1]];
	if ($c1+$c2<$min_depth)
	{
		$cache[$bam][$key] = "n/a";
		return "n/a";
	}
	
	//determine genotype
	$b1 = $keys[0];
	$b2 = ($c2>3 && $c2/($c1+$c2)>0.1) ? $keys[1] : $keys[0];
	
	$cache[$bam][$key] = "$b1/$b2";
	return "$b1/$b2";
}

//parse command line arguments
$parser = new ToolBase("kasp_ngs_match", "Tries to find matches for a KASP result in NGS datasets.");
$parser->addInfile("in", "KASP result file (_converted.tsv).", false);
$parser->addInfileArray("bams", "BAM files.", false);
$parser->addInt("min_depth", "Minimal depth for genotype determination in NGS data.", true, 10);
extract($parser->parse($argv));

//load KASP result
$data = [];
$file = file($in);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="") continue;
	
	if (starts_with($line, "#sample\t"))
	{
		$header = explode("\t", $line);
		continue;
	}
	
	$parts = explode("\t", $line);
	$sample = $parts[0];
	for ($i=1; $i<count($parts); ++$i)
	{
		$data[$sample][$header[$i]] = $parts[$i];
	}
}

//compare to bams
print "#kasp_sample\tbam\tsnps_kasp_and_ngs\tsnps_match_perc\n";
foreach($data as $sample => $var2geno)
{
	foreach($bams as $bam)
	{
		$c_kasp = 0;
		$c_both = 0;
		$c_match = 0;
		foreach($var2geno as $var => $g)
		{
			if ($g=="n/a") continue;
			++$c_kasp;
			
			list($rs, $chr, $pos, $ref, $alt) = explode("_", $var);
			$g2 = ngs_geno($bam, $chr, $pos, $ref, $min_depth);

			if ($g2=="n/a") continue;
			++$c_both;
					
			//normalize genotypes
			$g = explode("/", $g);
			sort($g);
			$g = implode("/", array_unique($g));
			$g_hom = (strlen($g)==1);
			
			$g2 = explode("/", $g2);
			sort($g2);
			$g2 = implode("/", array_unique($g2));
			$g2_hom = (strlen($g2)==1);
			
			if ($g==$g2) ++$c_match;
		}
		print "{$sample}\t{$bam}\t{$c_both}\t".number_format(100*$c_match/$c_both, 2)."\n";
	}
}

?>