<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

function determine_column($parts, $name)
{
	for($i=0; $i<count($parts);  ++$i)
	{
		if ($parts[$i]==$name) return $i;
	}
	
	trigger_error("Could not find column {$name}. In header columns:\n".implode("\n", $parts), E_USER_ERROR);
}

//init
$max_mb = $argv[1];
$classes = explode("/", $argv[2]);
$debug = false;

//parse input
$h = fopen2("php://stdin", "r");
while(!feof($h))
{
	$line = trim(fgets($h));
	if ($line=="") continue;
	
	$parts = explode("\t", $line);
	if (count($parts)<31) continue;
	
	//determine header indices
	if ($line[0]=="#")
	{
		$i_type = determine_column($parts, "Type");
		$i_name = determine_column($parts, "Name");
		$i_clnsig = determine_column($parts, "ClinicalSignificance");
		$i_pheno = determine_column($parts, "PhenotypeList");	
		$i_assembly = determine_column($parts, "Assembly");
		$i_chr = determine_column($parts, "Chromosome");
		$i_start = determine_column($parts, "Start");
		$i_end = determine_column($parts, "Stop");
		$i_id = determine_column($parts, "VariationID");		
		continue;
	}
	
	//check that it's a CNV
	$type = $parts[$i_type];
	if ($type=="copy number gain") 
	{
		$cn_type = "gain";
	}
	else if ($type=="copy number loss")
	{
		$cn_type = "loss";
	}
	else
	{
		if ($debug) print "Skipped line (no CNV): $line\n"; 
		continue;
	}
	
	//check that that GRCh38 coordinates are available
	$assembly = $parts[$i_assembly];
	$chr = null;
	if ($assembly=="GRCh38")
	{
		$chr = $parts[$i_chr];
		$start = $parts[$i_start];
		$end = $parts[$i_end];
	}
	$name = $parts[$i_name];
	if (starts_with($name, "GRCh38"))
	{
		if (preg_match("/.*\((.*):([0-9]+)-([0-9]+)\).*/", $name, $matches))
		{
			$chr = $matches[1];
			$start = $matches[2];
			$end = $matches[3];
		}
	}
	if (is_null($chr))
	{
		if ($debug) print "Skipped line (not GRCh38): $line\n"; 
		continue;
	}
	if ($start>=$end)
	{
		if ($debug) print "Skipped line (invalid range): $line\n"; 
		continue;
	}
	
	//check that it is annotated with a trend
	$clnsig = $parts[$i_clnsig];
	$class_hit = false;
	foreach($classes as $class)
	{
		if (contains($clnsig, $class)) $class_hit = true;
	}
	if (!$class_hit)
	{
		if ($debug) print "Skipped line (not pathogenic or benign): $line\n"; 
		continue;
	}
	
	//Skip variants bigger than max size
	if (($end-$start)/1000000.0 > $max_mb) continue;
	
	
	//extract exact CN state, if available
	list(, $cn) = explode(")x", $name.")xn/a");
	
	//output
	if (!starts_with($chr, "chr")) $chr = "chr".$chr;
	
	$id = $parts[$i_id];
	
	$phenotype = $parts[$i_pheno];
	if ($phenotype=="not provided") $phenotype = "n/a";

	print implode("\t", [ $chr, $start, $end, "$id [CN_TYPE={$cn_type} CN={$cn} CLNSIG=\"{$clnsig}\" PHENOTYPE=\"{$phenotype}\"]" ] )."\n";
}

	 
?>