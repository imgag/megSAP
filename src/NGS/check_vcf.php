<?php

/** 
	@page check_vcf
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("check_vcf", "Performs VCF file format check.");
$parser->addInfile("in",  "Input file in VCF format.", false);
//optional
$parser->addInt("limit", "Maximum number of lines that should be used to validate the vcf file (including header lines). '0' means all.", true, 10000);
extract($parser->parse($argv));

//unzip input VCF if necessary
$vcf_zipped = true;
if (ends_with($in, ".vcf"))
{
	$vcf_zipped = false;
}

//extract maximum number of lines
if ($limit<=0)
{
	$tmp_vcf = $in;
}
else
{
	if(!$vcf_zipped)
	{
		$tmp_vcf = $parser->tempFile(".vcf");
	}
	else
	{
		$tmp_vcf = $parser->tempFile(".vcf.gz");
	}
	$parser->exec("head", "$in -n $limit > $tmp_vcf", false);
}

//check that VCF file starts with format line (unless file is empty)
$handle = gzopen($tmp_vcf, "r");
$line = "";
while (!feof($handle) && $line=="")
{
	$line = trim(fgets($handle));
}
if ($line!="" && !starts_with($line, "##fileformat="))
{
	trigger_error("First line is not the VCF format line: $line", E_USER_ERROR);
}
gzclose($handle);

//check reference bases
$var_count = 0;
$handle = gzopen($tmp_vcf, "r");
while(!feof($handle))
{
	$line = trim(fgets($handle));
	if ($line=="" || $line[0]=="#") continue;

	list($chr, $start, , $ref, $alt) = explode("\t", $line);	
	$end = $start + strlen($ref) - 1;
	$ref2 = get_ref_seq($chr, $start, $end);
	if(strcasecmp($ref, $ref2)!=0)
	{
		trigger_error("Invalid reference sequence in variant $chr:$start $ref>$alt (ref should be: $ref2).", E_USER_ERROR);
	}
	
	++$var_count;
}
gzclose($handle);

//use SnpSIFT to check vcf-format - skip vcfcheck / vcf-output of freebayes not standard conform
//output of vcf contains .. for somatic variant files (?)
$stderr = "";
list(,$lines) = $parser->exec(get_path("SnpSift"), "vcfCheck $tmp_vcf",false);
foreach($lines as $line)
{
	$stderr .= str_replace("\t", " ",trim(nl_trim($line),". "));
}
if(!empty($stderr))	trigger_error("Invalid VCF-file. Output of vcfCheck: $stderr",E_USER_ERROR);

//check for valid MISO terms
$pattern = "ANN[*].EFFECT";
list($stdout, $stderr) = $parser->exec(get_path("SnpSift"), "extractFields $tmp_vcf '$pattern'", false);
$effects = explode("\t", str_replace("&", "\t", implode("\t", $stdout)));
$effects = array_map("trim", $effects);
$effects = array_unique($effects);
foreach($effects as $e)
{
	if ($e==$pattern || $e=="") continue;
	if(!Obo::isValidTermByName(repository_basedir()."/data/dbs/Ontologies/so-xp_3_0_0.obo",$e))
	{
		trigger_error("Invalid VCF-file. SnpEff effect '$e' is not valid Miso term in '$in'.", E_USER_ERROR);
	}
}

?>