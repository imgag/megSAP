<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//read gene id => names association from UCSC
$id2gene = array();
$handle = fopen("../UCSC/kgXref.txt", "r");
while(!feof($handle))
{
	$line = trim(fgets($handle));
	if ($line=="") continue;
	
	list ($id, , , , $name) = explode("\t", $line);
	$id2gene[$id] = strtoupper($name);
}
fclose($handle);

//load gene coordinates ranges from UCSC (with UTR and 20 flanking bases)
$gene2coord = array();
$handle = fopen("../UCSC/knownGene.txt", "r");
while(!feof($handle))
{
	$line = trim(fgets($handle));
	if ($line=="" || $line[0]=="#") continue;
	
	list($id, $chr, , $start, $end) = explode("\t", $line);
	if (!isset($id2gene[$id]))
	{
		trigger_error("Unknown UCSC identifier $id in knownGene.txt line: $line", E_USER_ERROR);
	}
	if ($start>=$end)
	{
		trigger_error("Invalid range $start-$end in knownGene.txt line: $line", E_USER_ERROR);
	}
	
	$tag = $chr."_".$id2gene[$id];
	$start -= 20;
	$end +=20;
	
	if (!isset($gene2coord[$tag]))
	{
		$gene2coord[$tag] = array($start, $end);
	}
	else
	{
		list($start_prev, $end_prev) = $gene2coord[$tag];
		$gene2coord[$tag] = array(min($start, $start_prev), max($end, $end_prev));
	}
}
fclose($handle);

//parse "mim2gene.txt" to get approved symbols
/*
# Mim Number	Type	Gene IDs	Approved Gene Symbols
100640	gene	216	ALDH1A1
100650	gene/phenotype	217	ALDH2
100660	gene	218	ALDH3A1
100670	gene	219	ALDH1B1
*/
$mim2gene = array();
$handle = fopen("mim2gene.txt", "r");
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="" || $line[0]=="#") continue;
	
	//only genes
	list($id, $type, , $gene) = explode("\t", $line);
	if (!contains($type, "gene")) continue;
	
	//only with approved symbols
	$gene = trim($gene);
	if ($gene=="") continue;
	
	$mim2gene[$id] = strtoupper($gene);
}
fclose($handle);

//parse "genemap.txt" for disorders
/*
# Sort	Month	Day	Year	Cyto Location	Gene Symbols	Confidence	Gene Name	MIM Number	Mapping Method	Comments	Phenotypes	Mouse Gene Symbol
1.1	5	13	13	1pter-p36.13	CTRCT8, CCV	P	Cataract, congenital, Volkmann type	115665	Fd	linked to Rh in Scottish family	Cataract 8, multiple types (2)	
1.2	9	25	01	1p36.23	ENO1, PPH, MPB1	C	Enolase-1, alpha	172430	S, F, R, REa		Enolase deficiency (1)	4(Eno1)
1.3	12	22	87	1pter-p36	ERPL1, HLM2	C	Endogenous retroviral pol gene-like sequence 1 (oncogene HLM2)	131190	REa, F			
1.4	4	14	11	1p36.11	HMGCL	P	3-hydroxy-3-methylglutaryl-Coenzyme A lyase	613898	REa, A		HMG-CoA lyase deficiency, 246450 (3)	4(Hmgcl)
1.5	4	30	15	1p36.33	AGRN, CMS8	P	Agrin	103320	REa		Myasthenic syndrome, congenital, 8, with pre- and postsynaptic defects, 615120 (3)	4(Agrn)
1.6	6	24	16	1p36.33	GNB1, MRD42	C	Guanine nucleotide-binding protein, beta polypeptide-1	139380	REa, A		Mental retardation, autosomal dominant 42, 616973 (3); Leukemia, acute lymphoblastic, somatic, 613065 (3)	4(Gnb1)
1.7	10	2	07	1p35.2	SDC3, SYND3, SDCN	C	Syndecan 3	186357	REc, R, H		{Obesity, association with}, 601665 (3)	4(Synd3)
1.8	8	28	98	1pter-p22.1	SAI1, MTS1, TFS1	C	Suppression of anchorage independence-1 (malignant transformation suppression-1)	154280	S, H			4(Tfs1)
1.9	12	6	16	1p36.33	ATAD3A, HAYOS	P	ATPase family, AAA domain-containing, member 3A	612316	REc	one family with AR inheritance reported	Harel-Yoon syndrome, 617183 (3)	4(Atad3)
*/
$handle = fopen("genemap.txt", "r");
print "#chr	start	end	omim\n";
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if($line=="" || $line[0]=="#") continue;
	
	$parts = explode("\t", $line);
	if (count($parts)<12)
	{
		print_r($parts);
		trigger_error("Line '$line' does not contain 12 tab-separated parts!", E_USER_ERROR);
	}
	
	//id
	$mim_id = trim($parts[8]);
	
	//extract genes (approved and not approved symbols - checked later)
	$genes = explode(",", $parts[5]);
	$genes = array_map("trim", $genes);
	if (isset($mim2gene[$mim_id])) $genes[] = $mim2gene[$mim_id];
	$genes = array_unique($genes);
	
	//confiddence not limbo
	$status = trim($parts[6]);
	if ($status=="C")
	{
		$status = "confirmed";
	}
	else if ($status=="P")
	{
		$status = "provisional";
	}
	else if ($status=="I")
	{
		$status = "inconsistent";
	}
	else continue;
	
	//disorders not empty
	$disorders = trim($parts[11]);
	$disorders = preg_replace("/, \d{6,6}/i", "", $disorders); 
	$disorders = trim(strtr($disorders, array("(1)"=>"", "(2)"=>"", "(3)"=>"", "(4)"=>"", "["=>"", "]"=>"", "{"=>"", "}"=>"", "?"=>"", ";"=>"|")));
	$disorders = trim(strtr($disorders, array(" , "=>",", " ,"=>",", ", "=>",")));
	$disorders = trim(strtr($disorders, array(" | "=>"|", " |"=>"|", "| "=>"|")));
	$disorders = trim(strtr($disorders, array(" "=>"_")));
	if ($disorders=="") continue;
	
	//chromosome
	$chr = "chr".trim($parts[0]);
	$chr = substr($chr, 0, strpos($chr, "."));
	if ($chr=="chr23") $chr="chrX";
	if ($chr=="chr24") $chr="chrY";
	
	//only with coordinates
	foreach($genes as $gene)
	{
		if(!isset($gene2coord[$chr."_".$gene])) continue;
		
		list($start, $end) = $gene2coord[$chr."_".$gene];
		print "$chr	$start	$end	{$mim_id}_[{$gene}_({$status})_{$disorders}]\n";
	}
}

fclose($handle);

?>