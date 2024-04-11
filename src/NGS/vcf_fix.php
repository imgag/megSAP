<?php 
/** 
	@page vcf_fix
	
	flags:
	--mosaic_mode: 
		- Doesn't remove variants with genotype '0/0' instead changes them to  het '0/1'.
	
	--longread_mode:
		- preserve phasing information
		- keeps different INFO/FORMAT values
	
	Fixes VCF file problems produced by Freebayes:
	- Merges duplicate heterozygous variants into one homozygous variant
	- Fixes 'AO' count for variants stemming from multi-allelic base variants (comma-separated list to single int)
	- Normalizes 'GT' to '0/1' or '1/1' (removes variants with genotype '0/0')
	- Removed variants with invalid REF characters ('N', ...) 
	- Removes unused fields in INFO/FORMAT columns
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//convert INFO data to associative array
function info_data($info)
{
	$info = explode(";", $info);
	$tmp = array();
	foreach($info as $entry)
	{
		if (!contains($entry, "="))
		{
			$tmp[$entry] = "";
		}
		else
		{
			list($key, $value) = explode("=", $entry, 2);
			$tmp[$key] = $value;
		}
	}
	
	return $tmp;
}

//convert FORMAT/SAMPLE data to associative array
function sample_data($format, $sample, $longread_mode=false)
{
	$tmp = array_combine(explode(":", $format), explode(":", $sample));
	
	if(!$longread_mode)
	{
		//normalize genotype string
		$tmp['GT'] = strtr($tmp['GT'], array("."=>"0", "|"=>"/"));
		if ($tmp['GT']=="1/0") $tmp['GT']="0/1";
	}
	
	return $tmp;
}

//write variant
function write($h_out, $var, $mosaic_mode=false, $longread_mode=false)
{
	// create format value column for each sample
	$format_header = "GT:DP:AO:GQ";
	if ($longread_mode) $format_header = "GT:DP:AF:GQ";
	$format_values = array();
	$all_wt = true;
	
	// Don't write variants with empty base info
	if ($var[0] == "" || $var[1] == "" || $var[2] == "" || $var[3] == "" || $var[4] == "")
	{
		return;
	}
	
	for ($i=9; $i < count($var); $i++) 
	{ 
		$sample= sample_data($var[8], $var[$i], $longread_mode);
		if($sample['GT']!="0/0" && $sample['GT']!=".|.")
		{
			//non wt/non-call variant -> keep in file
			$all_wt = false;
		}
		else if ($mosaic_mode)
		{
			// in mosaic mode keep wt variants, but change them into het variants.
			$sample['GT'] = "0/1";
			$all_wt = false;
		}
		
		// create output string
		if($longread_mode)
		{
			//set defaults: (some indels only contain GT column)
			if (!isset($sample['DP'])) $sample['DP'] = ".";
			if (!isset($sample['AF'])) $sample['AF'] = ".";
			if (!isset($sample['GQ'])) $sample['GQ'] = ".";

			$format_values[] = $sample['GT'].":".$sample['DP'].":".(($sample['AF'] == ".") ? "." :number_format($sample['AF'], 3)).":".(int)($sample['GQ']);
		}
		else
		{
			$format_values[] = $sample['GT'].":".$sample['DP'].":".(int)($sample['AO']).":".(int)($sample['GQ']);
		}
		
	}

	//skip wildtype variants
	if ($all_wt) 
	{
		return;
	}
	
	//skip variants with invalid reference characters
	if (!preg_match("/^[ACGT]+$/", $var[3]))
	{
		return;
	}
	
	//write base info
	fwrite($h_out, $var[0]."\t".$var[1]."\t".$var[2]."\t".$var[3]."\t".$var[4]."\t".number_format($var[5], 0, ".", "")."\t".$var[6]."\t");
	
	//write INFO
	$info = info_data($var[7]);

	if($longread_mode)
	{
		$info_output = array();
		if(isset($info["F"])) $info_output[] = "F";
		if(isset($info["P"])) $info_output[] = "P";
		// write '.' if INFO column is empty
		if(count($info_output) == 0) $info_output[] = ".";		fwrite($h_out, implode(";", $info_output)."\t");
	}
	else
	{
		fwrite($h_out, "MQM=".number_format($info['MQM'], 0, ".", "").";SAP=".number_format($info['SAP'], 0, ".", "").";SAR=".number_format($info['SAR'], 0, ".", "").";SAF=".number_format($info['SAF'], 0, ".", "").";ABP=".number_format($info['ABP'], 0, ".", "")."\t");
	}
	
	//write FORMAT/SAMPLE
	fwrite($h_out, $format_header."\t".implode("\t", $format_values)."\n");
}


//main loop

//CLI handling:
$mosaic_mode = false;
$longread_mode = false;
$format_header = "GT:DP:AO:GQ";
if ($argc != 1)
{
	if ($argv[1] == "--mosaic_mode")
	{
		$mosaic_mode = true;
	}
	elseif ($argv[1] == "--longread_mode")
	{
		$longread_mode = true;
		$format_header = "GT:DP:AF:GQ";		
	}
	else
	{
		fprintf(STDERR, "Unknown command line option(s) given '".$argv[1]."'.\nCurrently only the flags --mosaic_mode and --longread_mode are supported.");
		return 1;
	}
}


$ids = array("RO","GTI","NS","SRF","NUMALT","DP","QR","SRR","SRP","PRO","EPPR","DPB","PQR","ROOR","MQMR","ODDS","AN","RPPR","PAIREDR");
$var_last = array("","","","","","","","",$format_header,"0/0:0:0:0");
$h_in = fopen2("php://stdin", "r");
$h_out = fopen2("php://stdout", "w");
while(!feof($h_in))
{
	$line = trim(fgets($h_in));
	if ($line=="") continue;
	
	//write headers
	if (starts_with($line, "#"))
	{
		//fix header numbers
		if (starts_with($line, "##INFO=<ID="))
		{
			foreach($ids as $id)
			{
				if (starts_with($line, "##INFO=<ID=$id,"))
				{
					$line = str_replace("##INFO=<ID=$id,Number=1,","##INFO=<ID=$id,Number=.,", $line);
				}
			}
		}
		
		fwrite($h_out, $line."\n");
		continue;
	}
	
	//write content
	$var = explode("\t", $line); //chr, pos, id, ref, alt, qual, filter, info, format, sample_1, ..., sample_n
	
	//same variant twice => merge
	if ($var[0]==$var_last[0] && $var[1]==$var_last[1] && $var[3]==$var_last[3] && $var[4]==$var_last[4])
	{
		for ($i=9; $i < count($var); $i++) 
		{
			$sample = sample_data($var[8], $var[$i]);
			$sample_last = sample_data($var_last[8], $var_last[$i]);
			// check if merging can be performed and determine resulting genotype
			if (count($var) <= 10)
			{
				// single sample (only merging of 0/1 and 0/1 is allowed)
				if ($sample['GT']!="0/1" || $sample_last['GT']!="0/1")
				{
					//this happens sometimes when a large variant block overlaps two target region blocks (see https://github.com/ekg/freebayes/issues/351)
					trigger_error("Error: merging same variant {$var[0]}:{$var[1]} {$var[3]}>{$var[4]} with genotypes '".$sample['GT']."' and '".$sample_last['GT']." in single sample VCF is not possible!", E_USER_WARNING);
				}
				$gt = "1/1";
			}
			else
			{
				// multisample (also allow merging of 0/0 and 1/1 of the same sample)
				if (($sample['GT']=="0/1" && $sample_last['GT']=="1/1") || 
					(!$sample['GT']=="1/1" && $sample_last['GT']=="0/1") || 
					(!$sample['GT']=="1/1" && $sample_last['GT']=="1/1"))
				{
					trigger_error("Warning: merging same variant {$var[0]}:{$var[1]} {$var[3]}>{$var[4]} with genotypes '".$sample['GT']."' and '".$sample_last['GT']." in multisample VCF is not possible!", E_USER_WARNING);
				}

				if ($sample['GT']=="0/0" && $sample_last['GT']=="0/0")
				{
					// wt
					$gt = "0/0";
				}
				else if (($sample['GT']=="0/1" && $sample_last['GT']=="0/1") || ($sample['GT']=="1/1" && $sample_last['GT']=="0/0") || ($sample['GT']=="0/0" && $sample_last['GT']=="1/1"))
				{
					// hom
					$gt = "1/1";
				}
				else
				{
					// het
					$gt = "0/1";
				}
			}
			
			$dp = $sample['DP'];
			if ($longread_mode)
			{
				$af = $sample['AF'];
				if (contains($af, ",")) list($af) = explode(",", $af);
				$af_last = $sample_last['AF'];
				if (contains($af_last, ",")) list($af_last) = explode(",", $af_last);
				$af = min($dp, $af + $af_last);
			}
			else
			{
				$ao = $sample['AO'];
				if (contains($ao, ",")) list($ao) = explode(",", $ao);
				$ao_last = $sample_last['AO'];
				if (contains($ao_last, ",")) list($ao_last) = explode(",", $ao_last);
				$ao = min($dp, $ao + $ao_last);
			}
			$gq_last = $sample_last['GQ'];
			$gq = max($sample['GQ'], $gq_last);

			// determine genotype
			if ($longread_mode)
			{
				$var_last[$i] = "$gt:$dp:$af:$gq";
			}
			else
			{
				$var_last[$i] = "$gt:$dp:$ao:$gq";
			}
			
		} 
		// change FORMAT column at last to prevent FORMAT/SAMPLE key<->value errors 
		$var_last[8] = $format_header;
	}
	else
	{
		write($h_out, $var_last, $mosaic_mode, $longread_mode);
		$var_last = $var;
	}
}
write($h_out, $var_last, $mosaic_mode, $longread_mode);
fclose($h_in);
fclose($h_out);

?>