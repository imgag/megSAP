<?php
/** 
	@page validate_somatic_tumor_only
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("validate_somatic_tumor_only", "Validates the somatic variant calling performance on a tumor-only sample.");
$parser->addInfile("sample_var", "Germline variants before adding tumor variants (VCF.GZ). Make sure the variants are left-normalized.", false);
$parser->addInfile("tumor_var", "Germline variants of sample used to add artificial tumor variants (VCF.GZ). Make sure the variants are left-normalized.", false);
$parser->addInfile("roi", "Target region for the validation (high-confidence region of reference samples intersected with high-depth region of experiment.", false);
$parser->addInfileArray("calls", "Somatic variant call files (VCF.GZ). Sorted in order of ascending tumor content. Make sure the variants are left-normalized.", false);
$parser->addString("tum_content", "Comma seperated list for tumor content percentages for 'calls', e.g. '5,10,20,40'. If not given tum_content will be approximated by somatic variant AF. This may underestimate the actual tumor content.", true);
$parser->addOutfile("vars_details", "Output TSV file for variant details.", false);
$parser->addOutfile("af_details", "Output TSV file for AF-specific output.", false);
$parser->addEnum("caller", "The caller used to create the vcf files given in '-calls'", false, ["varscan", "deepsomatic"]);
$parser->addFlag("ignore_filters", "If set, variants with filter column entries are used. The default is to ignore variants with filter column entries.");
extract($parser->parse($argv));

//checks if the variant type matches
function type_matches($type, $tag)
{	
	if ($type=="SNVs" && is_indel($tag)) return false;
	if ($type=="InDels" && !is_indel($tag)) return false;
	
	return true;
}

function is_indel($tag)
{
	list($chr, $pos, $ref, $alt) = explode(" ", strtr($tag, ":>", "  "));
	return strlen($ref) > 1 || strlen($alt) > 1;
}

//load VCF variants. 'Caller' argument modifies information loaded for different callers (af and depth) and the reference files (genotype)
function load_vcf($filename, $roi, $caller)
{
	global $parser;
	global $ignore_filters;
	$c_filtered = 0;
	
	$output = array();
	list($lines) = $parser->execApptainer("htslib", "tabix", "--regions {$roi} {$filename}", [$roi], [dirname($filename)]);
	foreach($lines as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = explode("\t", $line."\t");
		$tag = "$chr:$pos $ref>$alt";
		
		//fix chr
		if (!starts_with($chr, "chr")) $chr = "chr".$chr;
		
		//parse gt/af
		if ($caller == "varscan")
		{
			list($dp, $af) = vcf_varscan2($format, $sample);

			//determine filter entries
			$filter = "";
			if (!$ignore_filters)
			{
				$filter = trim(strtr($filter, ["PASS"=>"", "."=>"", "freq-tum"=>"", ";"=>""]));
			}
			if ($filter!="") ++$c_filtered;
			
			$output[$tag] = array($af, $dp, $filter);
		}
		else if ($caller == "deepsomatic")
		{
			list($dp, $af) = vcf_deepvariant($format, $sample);
			//clear filter column
			$filter_replace = ["PASS"=>"", "freq-tum"=>"", ";"=>""];
			$filter = trim(strtr($filter, $filter_replace));
			if ($ignore_filters) $filter = "";
			
			$output[$tag] = array($af, $dp, $filter);
		}
		else if ($caller == "ref")
		{
			$gt = strtr(explode(":", $sample)[0], "/", "|");
			if ($gt=="1") $gt = "hom";
			else if ($gt=="1|1") $gt = "hom";
			else if ($gt=="1|0") $gt = "het";
			else if ($gt=="0|1") $gt = "het";
			else if ($gt==".|1") $gt = "het";
			else if ($gt=="1|.") $gt = "het";
			else  trigger_error("Unhandled genotype '{$gt}' in file '{$filename}'!", E_USER_ERROR);
			
			$output[$tag] = $gt;
		}
		else
		{
			trigger_error("Unknown 'caller' given in load_vcf(): $caller", E_USER_ERROR);
		}
	}
	return [$output, $c_filtered];
}

//create a sorted list with all unique variants in the given lists
function create_unique_variants_list($parser, $details_all)
{
	//create a unique variant list in VCF format
	$tags_done = [];
	$vcf = [];
	$vcf[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	foreach($details_all as $name => $file_vars)
	{
		foreach($file_vars as $tag => $result)
		{
			if (isset($tags_done[$tag])) continue;
			$tags_done[$tag] = true;
			
			list($chr, $pos, $ref, $alt) = explode(" ", strtr($tag, "-:>", "   "));
			$vcf[] = $chr."\t".$pos."\t.\t".$ref."\t".$alt."\t.\t.\t";
		}
	}
	$tmp = $parser->tempFile(".vcf");
	file_put_contents($tmp, implode("\n", $vcf));
	
	//sort variants
	$tmp2 = $parser->tempFile(".vcf");
	list($stdout) = $parser->execApptainer("ngs-bits", "VcfSort", "-in $tmp -out $tmp2");
	
	//create a tag list
	$output = [];
	foreach(file($tmp2) as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		list($chr, $pos, $id, $ref, $alt) = explode("\t", $line);
		$output[] = "$chr:$pos $ref>$alt";
	}
	
	return $output;
}

//***** MAIN SCRIPT *****\\

//verify given tumor content and parse expected heterogenic allele frequency
if ($tum_content != "")
{
	$parts = explode(",", $tum_content);
	if (count($parts) != count($calls)) trigger_error("There needs to be one tumor content per given calling file. ", E_USER_ERROR);
	
	$afs = [];
	foreach ($parts as $content)
	{
		$afs[] = (floatval($content) / 100) / 2; // base on expected AF of heterogenic variants
	}
}

//determine target region size
$roi_bases = bed_size($roi);
print "##ROI bases: {$roi_bases}\n";

//load germline variants
list($vars_germline, $filtered_germline) = load_vcf($sample_var, $roi, "ref");
print "##Variants of germline sample: ".count($vars_germline)." ({$filtered_germline} with filter entry)\n";

//load 'tumor' variants
list($vars_induced, $filtered_induced)  = load_vcf($tumor_var, $roi, "ref");
print "##Variants of 'tumor' sample: ".count($vars_induced)." ({$filtered_induced} with filter entry)\n";
foreach($vars_induced as $var => $gt)
{
	if (array_key_exists($var, $vars_germline))
	{
		unset($vars_induced[$var]);
	}
}
print "##  - after removing overlap with germline: ".count($vars_induced)."\n";

//benchmark
$details_all = []; //format: array($filename -> array($var -> [var_AF/"Missed"/"filtered"]))
print "#name\texpected\tTP\tTN\tFP\tFN\trecall/sensitivity\tprecision/ppv\tspecificity\n";
foreach($calls as $filename)
{
	//load calls
	$name = basename($filename, ".vcf.gz");
	list($vars, $filtered) = load_vcf($filename, $roi, $caller);	
	print "##Variant calls of '{$filename}': ".count($vars)." ({$filtered} with filter entry)\n";
	//remove expected germline variants for calls
	foreach($vars as $var => $gt)
	{
		if (array_key_exists($var, $vars_germline))
		{
			unset($vars[$var]);
		}
	}
	print "##  - after removing overlap with germline: ".count($vars)."\n";

	foreach(["", "SNVs", "InDels"] as $type)
	{
		$details = [];
		$c_expected = 0;
		$tp = 0;
		$fp = 0;
		$fn = 0;
		foreach($vars as $var => list($af, $dp, $filter))
		{
			if (!type_matches($type, $var)) continue;
			if (isset($vars_induced[$var]))
			{
				//flag filtered out
				if ($filter!="")
				{
					$af = "FILTERED";
				}
				else
				{
					++$tp;
				}
			}
			else if(!isset($vars_induced[$var]))
			{
				++$fp;
			}
			$details[$name][$var] = $af;
		}
		
		foreach($vars_induced as $var => $gt)
		{
			if (!type_matches($type, $var)) continue;
			++$c_expected;
			
			if (!isset($details[$name][$var]))
			{
				$details[$name][$var] = "MISSED";
				++$fn;
			}
		}
		
		$tn = $roi_bases - $c_expected - $fp;
		
		$recall = ($tp+$fn==0) ? "0.00000" : number_format($tp/($tp+$fn),5);
		$precision = ($tp+$fp==0) ? "0.00000" : number_format($tp/($tp+$fp),5);
		$spec = number_format($tn/($tn+$fp),5);
		
		print implode("\t", array(trim($name." ".$type), $c_expected, $tp, $tn, $fp, $fn, $recall, $precision, $spec))."\n";
		
		if ($type=="") $details_all[$name] = $details[$name];
	}
}


//create sorted and unique variant tag list
$vars_all = create_unique_variants_list($parser, $details_all);

//write header 
$var_detail_lines = [];
$header = "#variant\ttype";
foreach($details_all as $name => $tmp)
{
	$header .= "\t$name";
}
$var_detail_lines[] = $header;

//write variant details
foreach($vars_all as $var)
{
	if (isset($vars_germline[$var])) continue;
	
	$type = "ARTEFACT";
	if (isset($vars_induced[$var])) $type = "SOMATIC (".$vars_induced[$var].")";
	$line = "{$var}".(is_indel($var) ? " (INDEL)" : "")."\t{$type}";
	foreach($details_all as $name => $var_list)
	{
		$value = "";
		if (isset($var_list[$var]))
		{
			$value = $var_list[$var];
			if (is_numeric($value)) $value = number_format($value, 5);
		}
		$line .= "\t$value";
	}
	$var_detail_lines[] = $line;
}
file_put_contents($vars_details, implode("\n", $var_detail_lines));


if ($tum_content == "")
{
	$afs = array();
	//determine empirical AF of each column/experiment 
	foreach($var_detail_lines as $line)
	{
		if ($line[0]=="#") continue;

		$parts = explode("\t", $line);
		
		list($var, $type) = $parts;
		$is_indel = contains($var, "INDEL");
		if (!$is_indel && $type=="SOMATIC (het)")
		{
			for($col=2; $col<count($parts); ++$col)
			{
				$value = $parts[$col];
				if (is_numeric($value))
				{
					$afs[$col][] = $value;
				}
			}
		}
	}
	foreach($afs as $col => $values)
	{
		$afs[$col] = median($values);
	}
}


//calculate performance for given AF cutoff
$output_af = [];
foreach($afs as $i => $af)
{
	
	$name = array_keys($details_all)[$i];
	$output_af[] = "##tumor fraction {$name}: ".number_format(2.0*$af, 4);
}
$output_af[] = "#af\texpected\tTP\tFP\tFN\trecall/sensitivity\tprecision/ppv";
foreach([0.05, 0.075, 0.1, 0.125, 0.15, 0.2] as $min_af)
{
	foreach(["", "SNV", "INDEL"] as $curr_type)
	{
		$done = [];
		$tp = 0;
		$fp = 0;
		$fn = 0;
		foreach($var_detail_lines as $line)
		{
			if ($line[0]=="#") continue;
			$parts = explode("\t", $line);
			
			list($var, $type) = $parts;
			$is_indel = contains($var, "INDEL");
			if ($curr_type=="" || ($curr_type=="SNV" && !$is_indel) || ($curr_type=="INDEL" && $is_indel))
			{
				//count each variant only once
				if (isset($done[$var])) continue;
				
				//artefacts (FP)
				if ($type=="ARTEFACT")
				{
					for($col=2; $col<count($parts); ++$col)
					{
						$value = $parts[$col];
						if (is_numeric($value) && $value>=$min_af)
						{
							++$fp;
						}
					}
					
					$done[$var] = true;
					
					continue;
				}
				
				for($col=2; $col<count($parts); ++$col)
				{
					$af_expected = $afs[$col-2];
					if ($type=="SOMATIC (hom)") $af_expected *= 2.0;
					
					if ($af_expected>=$min_af)
					{
						$value = $parts[$col];
						if (!is_numeric($value)) //"MISSING", "FILTERED"
						{
							++$fn;
						}
						else
						{
							++$tp;
						}
						
						$done[$var] = true;
						
						break;
					}
				}
			}
		}

		$recall = ($tp+$fn==0) ? "n/a" : number_format($tp/($tp+$fn),5);
		$precision = ($tp+$fp==0) ? "n/a" : number_format($tp/($tp+$fp),5);
		$output_af[] = implode("\t", [ trim("{$min_af} {$curr_type}"), $tp+$fn, $tp, $fp, $fn, $recall, $precision ]);
	}
}
file_put_contents($af_details, implode("\n", $output_af));

?>
