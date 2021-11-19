<?php
/** 
	@page vcf_remove_duplicate_variants
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("vcf_remove_duplicate_variants", "Removes variants that are contained more than once.");
$parser->addInfile("in", "Input VCF file.", false);
$parser->addOutfile("out", "Output VCF file.", false);
extract($parser->parse($argv));


//determine duplicates
$duplicates = [];
$h = fopen2($in, "r");
$last_tag = "";
while(!gzeof($h))
{
	$line = fgets($h);
	if(trim($line)=="" || $line[0]=="#") continue;
	
	list($chr, $pos, $id, $ref, $alt) = explode("\t", $line);
	$tag = "$chr $pos $ref $alt";
	
	if ($tag==$last_tag)
	{
		$duplicates[$tag] = true;
	}
	
	$last_tag = $tag;
}
fclose($h);
print "Found ".count($duplicates)." duplicate variants.\n";

//remove duplicates
$ho = fopen2($out, "w");
$h = fopen2($in, "r");
while(!gzeof($h))
{
	$line = fgets($h);
	if(trim($line)=="" || $line[0]=="#")
	{
		fputs($ho, $line);
		continue;
	}
	
	list($chr, $pos, $id, $ref, $alt) = explode("\t", $line);
	$tag = "$chr $pos $ref $alt";
	
	if (isset($duplicates[$tag])) continue;
	
	fputs($ho, $line);
}
fclose($h);
fclose($ho);

?>