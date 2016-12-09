<?php
/** 
	@page vcf_check_leftaligned
*/

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";
require_once($basedir."Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("vcf_check_leftaligned", "\$Rev: 785 $", "Checks if variants in a VCF file are left-aligned.");
$parser->addInfile("in", "Input VCF file.", false);
//optional
$parser->addInt("lines", "The number of lines from the input VCF file to check", true, 100000);
extract($parser->parse($argv));

//(1) extract $lines lines from the head of the file
$head = $parser->tempFile(".vcf");
if (ends_with($in, ".gz"))
{
	exec2("zcat $in | head -n{$lines} | egrep -v \"CN[0-9]|INS:ME\" > $head");
}
else
{
	exec2("cat $in | head -n{$lines} | egrep -v \"CN[0-9]|INS:ME\" > $head");
}

//(2) left-align
$head_left = $parser->tempFile("_leftaligned.vcf");
exec2(get_path("ngs-bits")."VcfLeftNormalize -in $head -out $head_left", true);

//(3) compare
list($stdout) = exec2(get_path("ngs-bits")."SampleDiff -in1 $head -in2 $head_left -ei -window 0", false);
foreach($stdout as $line)
{
	print $line."\n";
}

?>

