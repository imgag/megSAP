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
$debug = false;

//parse input
$h = fopen2("php://stdin", "r");
while(!feof($h))
{
	$line = nl_trim(fgets($h));
	if ($line=="") continue;
	
	$parts = explode("\t", $line);
	if (count($parts)<17) continue;
	
	//determine header indices
	if ($line[0]=="#")
	{
		$i_gene = determine_column($parts, "#Gene Symbol");
		$i_loc = determine_column($parts, "Genomic Location");
		$i_class = determine_column($parts, "Haploinsufficiency Description");
		$i_pmid1 = determine_column($parts, "Haploinsufficiency PMID1");
		$i_pmid2 = determine_column($parts, "Haploinsufficiency PMID2");
		$i_pmid3 = determine_column($parts, "Haploinsufficiency PMID3");
		continue;
	}
	
	list($chr, $start, $end) = explode(":", strtr($parts[$i_loc], array("-"=>":", " "=>"")));
	$gene = trim($parts[$i_gene]);
	$class = trim($parts[$i_class]);
	if ($class=="No evidence available" || $class=="Dosage sensitivity unlikely") continue;
	$pmids = array($parts[$i_pmid1], $parts[$i_pmid2], $parts[$i_pmid3]);
	$pmids = array_map('trim', $pmids);
	$pmids = array_diff($pmids, array(""));

	print implode("\t", [ $chr, $start, $end, "$gene [HI={$class} PMIDS=".implode(",", $pmids)."]" ] )."\n";
}

	 
?>