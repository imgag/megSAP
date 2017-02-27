<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

/*
chrX	stdin	exon	2527306	2527522	.	+	.	gene_id "CD99P1"; transcript_id "NR_033380"; exon_number "1"; exon_id "NR_033380.1"; gene_name "CD99P1";
chrX	stdin	exon	2529037	2529079	.	+	.	gene_id "CD99P1"; transcript_id "NR_033380"; exon_number "2"; exon_id "NR_033380.2"; gene_name "CD99P1";
chrX	stdin	exon	2530191	2530259	.	+	.	gene_id "CD99P1"; transcript_id "NR_033380"; exon_number "3"; exon_id "NR_033380.3"; gene_name "CD99P1";
chrX	stdin	exon	2533857	2534010	.	+	.	gene_id "CD99P1"; transcript_id "NR_033380"; exon_number "4"; exon_id "NR_033380.4"; gene_name "CD99P1";
chrX	stdin	exon	2536781	2536825	.	+	.	gene_id "CD99P1"; transcript_id "NR_033380"; exon_number "5"; exon_id "NR_033380.5"; gene_name "CD99P1";
*/

$h1 = fopen("php://stdin", "r");
$h2 = fopen("php://stdout", "w");
while(!feof($h1))
{
	$line = trim(fgets($h1));
	if ($line=="") continue;
	
	$parts = explode("\t", $line); 
	
	//parse last column (annotations)
	$annos = array();
	foreach(explode(";", $parts[8]) as $anno)
	{
		$anno = trim($anno);
		if ($anno=="") continue;
		
		list($key, $value) = explode(" ", $anno);
		$annos[$key] = trim(substr($value, 1, -1));
	}
	
	//replace gene_id by transcript_id without suffix for gene duplications
	$trans = explode('_', $annos['transcript_id']);
	$annos['gene_id'] = $trans[0]."_".$trans[1];
	
	//replace last column
	$tmp = "";
	foreach($annos as $key => $value)
	{
		$tmp .= "$key \"$value\"; ";
	}
	$parts[8] = trim($tmp);
	
	//write output
	fwrite($h2, implode("\t", $parts)."\n");
}

fclose($h1);
fclose($h2);

?>