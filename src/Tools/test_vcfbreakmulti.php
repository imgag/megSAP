<?php
/** 
	@page test_vcfbreakmulti
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("test_vcfbreakmulti", "Tests vcflib 'vcfbreakmulti' function and ngs-bits 'VcfBreakMulti' function for a given vcf file.");
$parser->addInfile("in", "Input VCF file.", false);
$parser->addOutfile("out", "Output log.", false);
$parser->addString("tool", "Tool to use to break multi-allelic variants into several lines. Options are: ngs-bits, vcflib", false);

extract($parser->parse($argv));

//vcfbreakmulti output in tmp
$tmp_vcf = temp_file(".vcf");
$genome = genome_fasta("GRCh38");

//run vcflib vcfbreakmulti if vcflib is given as "tool"
if ($tool === "vcflib")
{
	$parser->execApptainer("vcflib", "vcfbreakmulti", "{$in} > {$tmp_vcf}", [$in]);
}
//run ngs-bits VcfBreakMulti if ngs-bits is given as "tool"
else if ($tool === "ngs-bits")
{
	$parser->execApptainer("ngs-bits", "VcfBreakMulti", "-in {$in} -out {$tmp_vcf}", [$in]);
}

//run ngs-bits VcfCheck on the output file and save the result
$parser->execApptainer("ngs-bits", "VcfCheck", "-in {$tmp_vcf} -out {$out} -ref {$genome}", [$genome], [dirname($out)]);

?>