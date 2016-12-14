<?php
/** 
	@page vcf_outliers
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("vcf_outliers", "Checks for a variant if there are outliers in the numeric info fields.");
$parser->addInfile("in", "Input VCF file.", false);
$parser->addString("var", "Variant position and base change, e.g. chr11:108121426:A:G", false);
extract($parser->parse($argv));

$stats = array();
$not_numeric = array();

$file = file($in);
$desc = array();
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="") continue;
	if ($line[0]=="#")
	{
		if (starts_with($line, "##INFO=<ID="))
		{
			$parts = explode(",", substr($line, 11, -1));
			$name = $parts[0];
			for($i=0; $i<count($parts); ++$i)
			{
				$part = $parts[$i];
				if (starts_with($part, "Description="))
				{
					$description = substr($part, 13);
					++$i;
					while($i<count($parts))
					{
						if (!contains($parts[$i], "="))
						{
							$description .= ",".$parts[$i];
							++$i;
							continue;
						}
						break;
					}
					break;
				}
			}
			$desc[$name] = strtr($description, array("\""=>""));
		}
		continue;
	}
	
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = explode("\t", $line);
	
	//check for variant
	if ("$chr:$pos:$ref:$alt"==$var)
	{
		$variant = $info;
		continue;
	}
	
	//accumulate numeric values
	$parts = explode(";", $info);
	foreach($parts as $part)
	{
		list($key, $value) = explode("=", $part);
		
		if (isset($not_numeric[$key])) continue;
		
		if (is_numeric($value))
		{
			$stats[$key][] = $value;
		}
		else
		{
			$not_numeric[$key] = true;
			if (isset($stats[$key]))
			{
				unset($stats[$key]);
			}
		}
		
	}
}

if (!isset($variant))
{
	trigger_error("Could not find variant '$var' in input file!", E_USER_ERROR);
}


//calculate statistics
print "#name\tvalue\tpercentile\tdescription\n";
$parts = explode(";", $variant);
foreach($parts as $part)
{
	list($key, $value) = explode("=", $part);
	$description = $desc[$key];
	
	if (isset($not_numeric[$key]))
	{
		print "$key\t$value\tnot numeric\t$description\n";
		continue;
	}
	
	//sort numeric
	$sort = $stats[$key];
	sort($sort,  SORT_NUMERIC);
	if (count($sort)<10)
	{
		print "$key\t$value\ttoo few dataponts\t$description\n";
		continue;
	}
	if ($sort[0]==$sort[count($sort)-1])
	{
		print "$key\t$value\tno variance\t$description\n";
		continue;
	}
	
	//determine percentile
	for ($i=0; $i<count($sort); ++$i)
	{
		if($value<=$sort[$i])
		{
			break;
		}
	}
	$sort = array_reverse($sort);
	for ($j=0; $j<count($sort); ++$j)
	{
		if($value>=$sort[$j])
		{
			break;
		}
	}
	$perc_start = number_format( 100*$i/count($sort),0);
	$perc_end = number_format(100-100*$j/count($sort),0);
	$perc = ($perc_start==$perc_end) ? "$perc_start" : "$perc_start-$perc_end";
	
	print "$key\t$value\t$perc\t$description\n";
}


?>

