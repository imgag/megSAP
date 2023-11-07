<?php
/** 
	@page diff_genome_builds
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("diff_genome_builds", "Compare variant lists of the same processed sample produced on GRCh37/GRCh38.");
$parser->addInfile("old", "GSvar file produced with GRCh37.", false);
$parser->addInfile("new", "GSvar file produced with GRCh38.", false);
extract($parser->parse($argv));

function load_gsvar($filename)
{
	$output = [];
	
	$headers = [];
	$file = file($filename);
	foreach($file as $line)
	{
		$line = nl_trim($line);
		if (trim($line)=="") continue;
		
		//header
		if ($line[0]=="#")
		{
			if (!starts_with($line, "##"))
			{
				$headers = explode("\t", $line);
			}
			continue;
		}
		
		//contents
		$parts = explode("\t", $line);
		$var = $parts[0].":".$parts[1]."-".$parts[2]." ".$parts[3].">".$parts[4];
		
		$annos = array_combine(array_slice($headers, 5), array_slice($parts, 5));
		$output[$var] = $annos;
	}
	
	return $output;
}

function max_num($str, $sep)
{
	$values = [];
	
	$parts = explode($sep, $str);
	foreach($parts as $part)
	{
		$part = trim($part);
		if (is_numeric($part)) $values[] = $part;
	}
	
	return count($values)==0 ? "" :max($values);
}

//load GSvar files
$old = load_gsvar($old);
print "##Variants old: ".count($old)."\n";
$new = load_gsvar($new);
print "##Variants new: ".count($new)."\n";

//match variants
$mapping = [];
$matches_dbsnp = 0;

$dbsnp2var = [];
foreach($new as $var => $annos)
{
	$dbsnp = trim($annos['dbSNP']);
	if ($dbsnp!="")
	{
		if (!isset($dbsnp2var[$dbsnp]))
		{
			$dbsnp2var[$dbsnp] = [];
		}
		
		$dbsnp2var[$dbsnp][] = $var;
	}
}

foreach($old as $var => $annos)
{	
	//dbSNP ID and base exchange, since dbSNP is not uniq, e.g. https://www.ncbi.nlm.nih.gov/snp/rs10547921
	$dbsnp = trim($annos['dbSNP']);
	if ($dbsnp!="")
	{
		if (isset($dbsnp2var[$dbsnp]))
		{
			foreach($dbsnp2var[$dbsnp] as $var2)
			{
				if (explode(" ", $var)[1]==explode(" ", $var2)[1])
				{
					$mapping[$var] = $var2;
					++$matches_dbsnp;
				}
			}
		}
	}
}
print "##Matching variants based on dbSNP: {$matches_dbsnp} (".number_format($matches_dbsnp * 100.0 / count($old), 2)."%)\n";

$text_cols = ["regulatory", "OMIM", "ClinVar", "HGMD", "RepeatMasker", "COSMIC", "MaxEntScan"]; //check text present/absent
$num_cols = ["1000g", "gnomAD", "phyloP", "CADD", "REVEL", "SpliceAI"]; //numeric comparison
$num_cols_sep = ["gnomAD_hom_hemi"=>",", "gnomAD_sub"=>","]; //numeric comparison (several values possible, max is taken)
$special_cols = ["quality"=>"depth", "coding_and_splicing"=>"impact", "NA12878_58"=>"genotype"]; //special handling

//init comparison result
print "#col\told\tnew_matching\tnew_matching_perc\tsimilarity\tcomment\n";
$cols_done = [];
$count_pos = [];
$count_match = [];
$similarity = [];
$comment = [];
foreach(array_merge($text_cols, $num_cols, array_keys($num_cols_sep), array_keys($special_cols)) as $col)
{
	$count_pos[$col] = 0;
	$count_match[$col] = 0;
	$similarity[$col] = "";
	$comment[$col] = "";
}

//compare text columns
foreach($text_cols as $col)
{
	$cols_done[$col] = true;
	
	foreach($mapping as $var => $var2)
	{
		$value = trim($old[$var][$col]);
		
		if ($value!="")
		{
			$count_pos[$col] += 1;
			
			$value2 = trim($new[$var2][$col]);
			$count_match[$col] += ($value2!="");
		}
	}
}

//compare numeric columns
foreach($num_cols as $col)
{
	$cols_done[$col] = true;
	
	$values = [];
	$values2 = [];
	foreach($mapping as $var => $var2)
	{
		$value = trim($old[$var][$col]);
		if ($value!="")
		{
			$count_pos[$col] += 1;
			
			$value2 = trim($new[$var2][$col]);
			if ($value2!="")
			{
				$count_match[$col] += 1;
			
				$values[] = $value;
				$values2[] = $value2;
			}
		}
	}
	$similarity[$col] = number_format(correlation($values, $values2), 2);
}

//compare numeric columns (with separator)
foreach($num_cols_sep as $col => $sep)
{
	$cols_done[$col] = true;
	
	$values = [];
	$values2 = [];
	foreach($mapping as $var => $var2)
	{
		$value = max_num($old[$var][$col], $sep);
		if ($value!="")
		{
			$count_pos[$col] += 1;
			
			$value2 = max_num($old[$var][$col], $sep);
			if ($value2!="")
			{
				$count_match[$col] += 1;
			
				$values[] = $value;
				$values2[] = $value2;
			}
		}
	}
	$similarity[$col] = number_format(correlation($values, $values2), 2);
}

//compare special columns
foreach($special_cols as $col => $com)
{
	$cols_done[$col] = true;
	$comment[$col] = $com;
	
	$values = [];
	$values2 = [];
	foreach($mapping as $var => $var2)
	{
		$value = trim($old[$var][$col]);
		$value2 = trim($new[$var2][$col]);
				
		if ($col=="quality")
		{
			$count_pos[$col] += 1;
			
			$entry = explode(";", $value)[1];
			$num = explode("=", $entry)[1];
			
			$entry2 = explode(";", $value2)[1];
			$num2 = explode("=", $entry2)[1];
			
			if ($num==$num2) $count_match[$col] += 1;
			
			$values[] = $num;
			$values2[] = $num2;
		}
		else if ($col=="coding_and_splicing")
		{
			$count_pos[$col] += 1;
			
			if (contains($value, ":HIGH:")) $num = 1.0;
			else if (contains($value, ":MODERATE:")) $num = 0.7;
			else if (contains($value, ":LOW:")) $num = 0.3;
			else $num = 0.0;
			
			if (contains($value2, ":HIGH:")) $num2 = 1.0;
			else if (contains($value2, ":MODERATE:")) $num2 = 0.7;
			else if (contains($value2, ":LOW:")) $num2 = 0.3;
			else $num2 = 0.0;
			
			if ($num==$num2) $count_match[$col] += 1;
			
			$values[] = $num;
			$values2[] = $num2;
		
		}
		else if ($col=="NA12878_58")
		{
			$count_pos[$col] += 1;
			$count_match[$col] += ($value==$value2);
		}
	}
	
	if ($col=="quality" || $col=="coding_and_splicing")
	{	
		$similarity[$col] = number_format(correlation($values, $values2), 2);
	}
}

//output
foreach($count_pos as $col => $c_pos)
{
	$c_match = $count_match[$col];
	print "{$col}\t{$c_pos}\t{$c_match}\t".number_format(100.0 * $c_match / $c_pos, 2)."\t".$similarity[$col]."\t".$comment[$col]."\n";
}

//unchecked columns
$columns = array_keys($old[$var]);
foreach($columns as $col)
{
	if (!isset($cols_done[$col]))
	{
		print "##Unchecked column: $col\n";
	}
}

?>