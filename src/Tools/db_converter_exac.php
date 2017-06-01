<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$subpopulations = array("AFR"=>true, "AMR"=>true, "EAS"=>true, "NFE"=>true, "SAS"=>true);

$handle = fopen("php://stdin", "r");
while($line = fgets($handle))
{
	$line = trim($line);
	if ($line=="") continue;
	if($line[0]=="#")
	{
		print $line."\n";
		continue;
	}
	
	$parts = explode("\t", $line);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = $parts;
	$info = explode(";", $info);
	
	$info_new = array();
	$ac = array();
	$an = array();
	foreach ($info as $entry)
	{
		//keep
		if (starts_with($entry, "AC_Hom=") || starts_with($entry, "Hom_NFE=") || starts_with($entry, "Hom_AFR="))
		{
			$info_new[] = $entry;
		}

		//calcualte overall AF from adjusted AC/AN
		else if (starts_with($entry, "AC_Adj="))
		{
			$ac['ALL'] = substr($entry, 7);
		}
		else if (starts_with($entry, "AN_Adj="))
		{
			$an['ALL'] = substr($entry, 7);
		}
		
		//subpopulation AF
		else if (starts_with($entry, "AC_") && $entry[6]=="=")
		{
			$pop = substr($entry, 3, 3);
			if (!isset($subpopulations[$pop])) continue;
			
			$ac[$pop] = substr($entry, 7);
		}
		else if (starts_with($entry, "AN_") && $entry[6]=="=")
		{
			$pop = substr($entry, 3, 3);
			if (!isset($subpopulations[$pop])) continue;
			
			$an[$pop] = substr($entry, 7);
		}
	}
	
	if (count($ac)!=count($an))
	{
		print $line;
		print_r($ac);
		print_r($an);
		trigger_error("Allele count and allele number array have different sizes!", E_USER_ERROR);
	}
	
	//print $line;
	foreach($ac as $key => $value)
	{
		if (!isset($an[$key]))
		{
			print $line;
			print_r($ac);
			print_r($an);
			trigger_error("Key '$key' not set in allle number array!", E_USER_ERROR);
		}
		
		if ($key=='ALL') //overall
		{
			$info_new[] = 'AF='.($an[$key]<200 ? '0.0' : number_format($value/$an[$key], 5));
		}
		else //subpopulations
		{
			$info_new[] = "AF_{$key}=".($an[$key]<200 ? "0.0" : number_format($value/$an[$key], 5));
		}
	}
	print "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t".implode(";", $info_new)."\n";
}
fclose($handle);

?>