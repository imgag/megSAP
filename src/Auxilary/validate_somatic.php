<?php
/** 
	@page validate_somatic
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("validate_somatic", "Validates the somatic variant calling performance.");
$parser->addInfile("normal", "Expected germline variants (VCF.GZ).", false);
$parser->addInfile("tumor", "Expected tumor variants  (VCF.GZ).", false);
$parser->addInfile("roi", "Target region for the validation (high-confidence region of reference samples intersected with high-depth region of experiment.", false);
$parser->addInfileArray("calls", "Somatic variant call files (VCF.GZ). Sorted in order of ascending tumor content", false);
$parser->addString("tum_content", "Tumor content in the samples used for the calls given as comma seperated string '5,10,20,40'. If not given tum_content will be approximated from variant af. This may underestimate actual content", true);
$parser->addOutfile("vars_details", "Output TSV file for variant details.", false);
$parser->addOutfile("af_details", "Output TSV file for AF-specific output.", false);
$parser->addEnum("caller", "The caller used to create the vcf files given in '-calls'", false, ["strelka2", "dragen", "deepsomatic"]);
$parser->addFlag("with_low_evs", "Ignores 'LowEVS' filter column entry for strelka2 callings");
$parser->addFlag("ignore_filters", "Ignores filter column entries.");
extract($parser->parse($argv));

//returns if the variant is an InDels
function is_indel($tag)
{
	list($chr, $pos, $ref, $alt) = explode(" ", strtr($tag, ":>", "  "));
	return strlen($ref) > 1 || strlen($alt) > 1;
}

//load VCF variants. 'Caller' argument can modify information loaded for different callers (af and depth) and the reference files (genotype)
function load_vcf($filename, $roi, $caller)
{
	global $with_low_evs;
	global $parser;
	global $ignore_filters;
	
	$output = array();
	list($lines) = $parser->execApptainer("htslib", "tabix", "--regions {$roi} {$filename}", [$roi], [dirname($filename)]);
	foreach($lines as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		if ($caller == "deepsomatic") list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample_tumor) = explode("\t", $line."\t");
		else list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample_normal, $sample_tumor) = explode("\t", $line."\t");
		//fix chr
		if (!starts_with($chr, "chr")) $chr = "chr".$chr;
		
		//fix gt/af
		$tag = "$chr:$pos $ref>$alt";
		
		if ($caller == "strelka2")
		{
			list($tumor_dp, $tumor_af) = is_indel($tag) ? vcf_strelka_indel($format, $sample_tumor) : vcf_strelka_snv($format, $sample_tumor, $alt);
			
			//clear filter column
			$filter_replace = ["PASS"=>"", "freq-tum"=>"", ";"=>""];
			if ($with_low_evs) $filter_replace ["LowEVS"] = "";
			$filter = trim(strtr($filter, $filter_replace));
			if ($ignore_filters) $filter = "";
			
			$output[$tag] = array($tumor_af, $tumor_dp, $filter);
		}
		else if ($caller == "dragen")
		{
			list($tumor_dp, $tumor_af) = vcf_dragen_var($format, $sample_tumor);
			
			//clear filter column
			$filter_replace = ["PASS"=>"", "."=>"", "freq-tum"=>"", ";"=>""];
			$filter = trim(strtr($filter, $filter_replace));
			if ($ignore_filters) $filter = "";
			
			$output[$tag] = array($tumor_af, $tumor_dp, $filter);
		}
		else if ($caller == "deepsomatic")
		{
			list($tumor_dp, $tumor_af) = vcf_deepvariant($format, $sample_tumor);

			//clear filter column
			$filter_replace = ["PASS"=>"", "freq-tum"=>"", ";"=>""];
			$filter = trim(strtr($filter, $filter_replace));
			if ($ignore_filters) $filter = "";
			
			$output[$tag] = array($tumor_af, $tumor_dp, $filter);
		}
		else if ($caller == "ref")
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
		else
		{
			trigger_error("Unknown 'caller' given in load_vcf(): $caller", E_USER_ERROR);
		}
	}
	return $output;
}

//checks if the variant type matches
function type_matches($type, $tag)
{	
	if ($type=="SNVs" && is_indel($tag)) return false;
	if ($type=="InDels" && !is_indel($tag)) return false;
	
	return true;
}

//create a sorted list with all unique variants in the given lists
function create_unique_variants_list($parser, $details_all)
{
	$vars_all = [];
	$vcf = [];
	$vcf[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	foreach($details_all as $name => $file_vars)
	{
		foreach($file_vars as $var => $result)
		{
			list($chr, $pos, $ref, $alt) = explode(" ", strtr($var, "-:>", "   "));
			$vcf[] = $chr."\t".$pos."\t.\t".$ref."\t".$alt."\t100\tPASS\t";
		}
	}
	$tmp = $parser->tempFile(".vcf");
	file_put_contents($tmp, implode("\n", $vcf));
	$tmp2 = $parser->tempFile(".vcf");
	list($stdout) = $parser->execApptainer("ngs-bits", "VcfSort", "-in $tmp -out $tmp2 && sort --uniq $tmp2");
	foreach($stdout as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		list($chr, $pos, $id, $ref, $alt) = explode("\t", $line);
		$vars_all[] = "$chr:$pos $ref>$alt";
	}
	
	return $vars_all;
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
$vars_germline = load_vcf($normal, $roi, "ref");
print "##Germline variants: ".count($vars_germline)."\n";

//load expected somatic variants
$vars_somatic = load_vcf($tumor, $roi, "ref");
print "##Tumor variants: ".count($vars_somatic)."\n";
foreach($vars_germline as $var => $gt)
{
	if (array_key_exists($var, $vars_somatic))
	{
		unset($vars_somatic[$var]);
	}
}
print "##Tumor variants after removing overlap with germline: ".count($vars_somatic)."\n";

//benchmark
$details_all = []; //format: array($filename -> array($var -> [var_AF/"Missed"/"filtered"]))
print "#name\texpected\tTP\tTN\tFP\tFN\trecall/sensitivity\tprecision/ppv\tspecificity\n";
foreach($calls as $filename)
{
	$name = basename($filename, "_var.vcf.gz");
	$vars = load_vcf($filename, $roi, $caller);
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
			if (isset($vars_somatic[$var]))
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
			else if(!isset($vars_germline[$var]))
			{
				++$fp;
			}
			$details[$name][$var] = $af;
		}
		
		foreach($vars_somatic as $var => $gt)
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


//create sorted variant list
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
	if (isset($vars_somatic[$var])) $type = "SOMATIC (".$vars_somatic[$var].")";
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

		$recall = ($tp+$fn==0) ? "n/a" : number_format($tp/($tp+$fn) * 100,2);
		$precision = ($tp+$fp==0) ? "n/a" : number_format($tp/($tp+$fp) * 100,2);
		$output_af[] = implode("\t", [ trim("{$min_af}% {$curr_type}"), $tp+$fn, $tp, $fp, $fn, $recall."%", $precision."%" ]);
	}
}
file_put_contents($af_details, implode("\n", $output_af));

?>
