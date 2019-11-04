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
	//skip wildtype variants
	$sample = sample_data($var[8], $var[9]);
	if ($sample['GT']=="0/0") return;
	
	//write base info
	fwrite($h_out, $var[0]."\t".$var[1]."\t".$var[2]."\t".$var[3]."\t".$var[4]."\t".number_format($var[5], 0, ".", "")."\t".$var[6]."\t");
	
	//write INFO
	$info = info_data($var[7]);
	fwrite($h_out, "MQM=".number_format($info['MQM'], 0, ".", "")."\t");

	//write FORMAT/SAMPLE
	fwrite($h_out, "GT:DP:AO\t".$sample['GT'].":".$sample['DP'].":".(int)($sample['AO'])."\n");
}


//main loop
$ids = array("RO","GTI","NS","SRF","NUMALT","DP","QR","SRR","SRP","PRO","EPPR","DPB","PQR","ROOR","MQMR","ODDS","AN","RPPR","PAIREDR");
$var_last = array("","","","","","","","","GT","0/0");
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
	$var = explode("\t", $line); //chr, pos, id, ref, alt, qual, filter, info, format, sample
	
	//do not convert multi-sample files
	if (count($var)>10)
	{
		fwrite($h_out, "$line\n");
		continue;
	}
	
	//same variant twice => merge
	if ($var[0]==$var_last[0] && $var[1]==$var_last[1] && $var[3]==$var_last[3] && $var[4]==$var_last[4])
	{
		$sample = sample_data($var[8], $var[9]);
		$sample_last = sample_data($var[8], $var[9]);
		if ($sample['GT']!="0/1" || $sample_last['GT']!="0/1")
		{
			//this happens sometimes when a large variant block overlaps two target region blocks (see https://github.com/ekg/freebayes/issues/351)
			$stderr = fopen2('php://stderr', 'w');
			fwrite($stderr, "Warning: merging same variant {$var[0]}:{$var[1]} {$var[3]}>{$var[4]} with genotypes '".$sample['GT']."' and '".$sample_last['GT']."'\n");
			fclose($stderr);
		}
		$dp = $sample['DP'];
		$ao = $sample['AO'];
		if (contains($ao, ",")) list($ao) = explode(",", $ao);
		$ao_last = $sample_last['AO'];
		if (contains($ao_last, ",")) list($ao_last) = explode(",", $ao_last);
		$ao = min($dp, $ao + $ao_last);
		$var_last[8] = "GT:DP:AO";
		$var_last[9] = "1/1:$dp:$ao";
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