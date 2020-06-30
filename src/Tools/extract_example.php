<?php
/** 
	@page extract_example 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("extract_example", "Extracts a region from a sample e.g to create demo data.");
$parser->addInfile("in", "BAM file of the sample.", false);
$parser->addString("out", "Output name.", false);
$parser->addString("reg", "Region to extract in format chr:start-end.", false);
extract($parser->parse($argv));

//init
$samtools = get_path("samtools");

//make sure output folder exists
$folder = $out."/";
exec2("mkdir -p $folder");

//store region as BED file
$roi = $parser->tempFile(".bed");
$reg = strtr($reg, array(","=>""));
if(!starts_with($reg, "chr")) $reg = "chr".$reg;
file_put_contents($roi, strtr($reg, array(":"=>"\t", "-"=>"\t")));

//extract reads from BAM
$out_bam = "{$folder}/{$out}.bam";
exec2("{$samtools} view -h -L {$roi} -M {$in} | {$samtools} view -Sb > {$out_bam}");
$parser->indexBam($out_bam, 4);

//write IGV session file
$session = array();
$session[] = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>";
$session[] = "<Session genome=\"1kg_v37\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"{$reg}\" path=\"igv_session.xml\" version=\"8\">";
$session[] = "    <Resources>";
$session[] = "        <Resource path=\"{$out}.bam\"/>";
$session[] = "    </Resources>";
$session[] = "    <Panel name=\"Panel1507276496064\">";
$session[] = "        <Track displayMode=\"COLLAPSED\" id=\"{$out}.bam_coverage\" name=\"{$out}.bam Coverage\" showReference=\"false\" snpThreshold=\"0.05\">";
$session[] = "        </Track>";
$session[] = "        <Track displayMode=\"EXPANDED\" id=\"{$out}.bam\" name=\"{$out}.bam\">";
$session[] = "            <RenderOptions viewPairs=\"true\"/>";
$session[] = "        </Track>";
$session[] = "    </Panel>";
$session[] = "    <PanelLayout dividerFractions=\"0.01,0.9\"/>";
$session[] = "</Session>";
$out_igv = "{$folder}/igv_session.xml";
file_put_contents($out_igv, implode("\n", $session));


?>