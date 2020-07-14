<?php
/** 
	@page validate_somatic
	//TODO: split in SNV/Indel
	//TODO: Add min_dp parameter
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);


//parse command line arguments
$parser = new ToolBase("validate_somatic", "Validates the somatic variant calling performance.");
$parser->addInfile("normal", "Expected germline variants (VCF.GZ).", false);
$parser->addInfile("tumor", "Expected tumor variants  (VCF.GZ).", false);
$parser->addInfile("roi", "Target region for the validation.", false);
$parser->addInfileArray("calls", "Somatic variant call files (VCF.GZ).", false);
$parser->addOutfile("vars_details", "Output TSV file for variant details.", false);
extract($parser->parse($argv));


function load_vcf($filename, $roi, $strelka)
{
	$output = array();
	list($lines) = exec2("tabix --regions {$roi} {$filename}");
	foreach($lines as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample_normal, $sample_tumor) = explode("\t", $line."\t");
		$is_indel = strlen($ref) > 1 || strlen($alt) > 1;
		
		//fix chr
		if (!starts_with($chr, "chr")) $chr = "chr".$chr;
		
		//fix gt/af
		$tag = "$chr:$pos $ref>$alt";
		if ($strelka)
		{
			list($tumor_dp, $tumor_af) = $is_indel ? vcf_strelka_indel($format, $sample_tumor) : vcf_strelka_snv($format, $sample_tumor, $alt);
			$output[$tag] = $tumor_af;
		}
		else
		{
			$gt = strtr(explode(":", $sample_normal)[0], "/", "|");
			if ($gt=="1") $gt = "hom";
			else if ($gt=="1|1") $gt = "hom";
			else if ($gt=="1|0") $gt = "het";
			else if ($gt=="0|1") $gt = "het";
			else if ($gt==".|1") $gt = "het";
			else if ($gt=="1|.") $gt = "het";
			else  trigger_error("Unhandled genotype '{$gt}' in file '{$filename}'!", E_USER_ERROR);
			
			$output[$tag] = $gt;
		}
	}
	return $output;
}

//load germline variants
$vars_germline = load_vcf($normal, $roi, false);
print "##Germline variants: ".count($vars_germline)."\n";

//load somatic variants
$tmp = $parser->tempFile(".bed");
exec2("cat {$roi} | tr -d 'chr' > {$tmp}");
$vars_somatic = load_vcf($tumor, $tmp, false);
list($stdout, $stderr) = exec2("BedInfo -in {$tmp} | grep Bases");
$roi_bases = trim(explode(":", $stdout[0])[1]);
print "##ROI bases: {$roi_bases}\n";

print "##Tumor variants: ".count($vars_somatic)."\n";
foreach($vars_germline as $var => $gt)
{
	if (isset($vars_somatic[$var]))
	{
		unset($vars_somatic[$var]);
	}
}
print "##Tumor variants after removing overlap with germline: ".count($vars_somatic)."\n";

//benchmark
$details = [];
print "#name\texpected\tTP\tTN\tFP\tFN\trecall/sensitivity\tprecision/ppv\tspecificity\n";
foreach($calls as $filename)
{
	$name = basename($filename, ".vcf.gz");
	$vars = load_vcf($filename, $roi, true);
	
	$tp = 0;
	$fp = 0;
	$fn = 0;
	foreach($vars as $var => $af)
	{
		if (isset($vars_somatic[$var]))
		{
			++$tp;
		}
		else if(!isset($vars_germline[$var]))
		{
			++$fp;
		}
		$details[$name][$var] = $af;
	}
	
	foreach($vars_somatic as $var => $gt)
	{
		if (!isset($details[$name][$var]))
		{
			$details[$name][$var] = "MISSED";
			++$fn;
		}
	}
	
	$c_expected = count($vars_somatic);
	$tn = $roi_bases - $c_expected - $fp;
	
	$recall = number_format($tp/($tp+$fn),5);
	$precision = number_format($tp/($tp+$fp),5);
	$spec = number_format($tn/($tn+$fp),5);
	
	print implode("\t", array(basename($filename, ".vcf.gz"), $c_expected, $tp, $tn, $fp, $fn, $recall, $precision, $spec))."\n";
}

//create sorted variant list
$vars_all = [];
$vcf = [];
$vcf[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
foreach($details as $name => $tmp)
{
	foreach($tmp as $var => $result)
	{
		list($chr, $pos, $ref, $alt) = explode(" ", strtr($var, "-:>", "   "));
		$vcf[] = $chr."\t".$pos."\t.\t".$ref."\t".$alt."\t100\tPASS\t";
	}
}
$tmp = $parser->tempFile(".vcf");
file_put_contents($tmp, implode("\n", $vcf));
$tmp2 = $parser->tempFile(".vcf");
list($stdout) = exec2(get_path("ngs-bits")."VcfSort -in $tmp -out $tmp2 && sort --uniq $tmp2", true);
foreach($stdout as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	list($chr, $pos, $id, $ref, $alt) = explode("\t", $line);
	$vars_all[] = "$chr:$pos $ref>$alt";
}

//write details
$output = [];
$header = "#variant\ttype";
foreach($details as $name => $tmp)
{
	$header .= "\t$name";
}
$output[] = $header;
foreach($vars_all as $var)
{
	if (isset($vars_germline[$var])) continue;
	$type = "ARTEFACT";
	if (isset($vars_somatic[$var])) $type = "SOMATIC (".$vars_somatic[$var].")";
	$line = "{$var}\t{$type}";
	foreach($details as $name => $tmp)
	{
		$value = "";
		if (isset($tmp[$var])) $value = $tmp[$var];
		$line .= "\t$value";
	}
	$output[] = $line;
}

file_put_contents($vars_details, implode("\n", $output));

?>
