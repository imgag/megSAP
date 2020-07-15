<?php 
/** 
	@page vcf_fix
	
	Fixes VCF file problems produced by Freebayes:
	- Merges duplicate heterozygous variants into one homozygous variant
	- Fixes 'AO' count for variants stemming from multi-allelic base variants (comma-separated list to single int)
	- Normalizes 'GT' to '0/1' or '1/1' (removes variants with genotype '0/0')
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
function sample_data($format, $sample)
{
	$tmp = array_combine(explode(":", $format), explode(":", $sample));
	
	//normalize genotype string
	$tmp['GT'] = strtr($tmp['GT'], array("."=>"0", "|"=>"/"));
	if ($tmp['GT']=="1/0") $tmp['GT']="0/1";
	
	return $tmp;
}

//write variant
function write($h_out, $var)
{
	// create format value column for each sample
	$format_values = array();
	$all_wt = true;
	for ($i=9; $i < count($var); $i++) 
	{ 
		$sample= sample_data($var[8], $var[$i]);
		if($sample['GT']!="0/0")
		{
			//non wt variant -> keep in file
			$all_wt = false;
		}
		// create output string
		$format_values[] = $sample['GT'].":".$sample['DP'].":".(int)($sample['AO']);
	}

	//skip wildtype variants
	if ($all_wt) return;
	
	//write base info
	fwrite($h_out, $var[0]."\t".$var[1]."\t".$var[2]."\t".$var[3]."\t".$var[4]."\t".number_format($var[5], 0, ".", "")."\t".$var[6]."\t");
	
	//write INFO
	$info = info_data($var[7]);
	fwrite($h_out, "MQM=".number_format($info['MQM'], 0, ".", "").";SAP=".number_format($info['SAP'], 0, ".", "").";ABP=".number_format($info['ABP'], 0, ".", "")."\t");

	//write FORMAT/SAMPLE
	fwrite($h_out, "GT:DP:AO\t".implode("\t", $format_values)."\n");
}


//main loop
$ids = array("RO","GTI","NS","SRF","NUMALT","DP","QR","SRR","SRP","PRO","EPPR","DPB","PQR","ROOR","MQMR","ODDS","AN","RPPR","PAIREDR");
$var_last = array("","","","","","","","","GT:DP:AO","0/0:0:0");
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
			$ao = $sample['AO'];
			if (contains($ao, ",")) list($ao) = explode(",", $ao);
			$ao_last = $sample_last['AO'];
			if (contains($ao_last, ",")) list($ao_last) = explode(",", $ao_last);
			$ao = min($dp, $ao + $ao_last);

			// determine genotype
			$var_last[$i] = "$gt:$dp:$ao";
		} 
		// change FORMAT column at last to prevent FORMAT/SAMPLE key<->value errors 
		$var_last[8] = "GT:DP:AO";
	}
	else
	{
		write($h_out, $var_last);
		$var_last = $var;
	}
}
write($h_out, $var_last);
fclose($h_in);
fclose($h_out);

?>