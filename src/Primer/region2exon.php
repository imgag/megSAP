<?php
/** 
	@page region2exon
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("region2exon", "Looks up exons for given regions.");
$parser->addString("reg",  "Input region in the format 'chr1:151137550-151138268'.", false);
$parser->addOutfile("out",  "Output BED file.", false);

extract($parser->parse($argv));

$output = array();

//extract matching exons
list($c, $s, $e) = preg_split("/[\s:-]+/", $reg);
$s -= 1; //subtract one from start to create correct BED format
$pipeline = [
	["echo", "'$c\t$s\t$e'"],
	[get_path("ngs-bits")."BedIntersect", "-in2 ".get_path("GRCh37_data_folder")."/dbs/UCSC/exons.bed -mode in2" ]
];
list($exons) = $parser->execPipeline($pipeline, "extrac matching exons");
$hits = array();
foreach ($exons as $exon)
{
	$exon = trim($exon);
	if ($exon=="") continue;
	$hits[] = array_slice(explode("\t", $exon), 1); //start, end, gene_exon, strand
}

//no matches
if (count($hits)==0)
{
	$output[] = "Not within an exon";
}

//mutiple matches: check whether there are mutiple distinct exons
for ($i=1; $i<count($hits); ++$i)
{
	if ($hits[$i][0]!=$hits[0][0] || $hits[$i][1]!=$hits[0][1])
	{
		$output[] = "Non-synonymous exons in region";
		break;
	}
}

//write all exons
foreach ($hits as $hit)
{
	$output[] = $c."\t".implode("\t", $hit);
}

file_put_contents($out, implode("\n", $output));
?>