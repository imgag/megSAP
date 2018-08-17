<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

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
	
	//extract genes
	$genes = explode(",", $parts[5]);
	$genes = array_map("trim", $genes);
	if (isset($mim2gene[$mim_id])) $genes[] = $mim2gene[$mim_id];
	$genes = array_unique($genes);
	
	//confidence not limbo
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
	$disorders = trim(strtr($disorders, array("(1)"=>"", "(2)"=>"", "(3)"=>"", "(4)"=>"", "["=>"", "]"=>"", "{"=>"", "}"=>"", "?"=>"", ";"=>"|")));
	$disorders = trim(strtr($disorders, array(" , "=>",", " ,"=>",", ", "=>",")));
	$disorders = trim(strtr($disorders, array(" | "=>"/", " |"=>"/", "| "=>"/")));
	$disorders = trim(strtr($disorders, array(" "=>"_")));
	if ($disorders=="") continue;
	
	//chromosome
	$chr_omim = "chr".trim($parts[0]);
	$chr_omim = substr($chr_omim, 0, strpos($chr_omim, "."));
	if ($chr_omim=="chr23") $chr_omim="chrX";
	if ($chr_omim=="chr24") $chr_omim="chrY";
	
	//convert genes to approved symbols
	$genes_approved = array();
	list($stdout) = exec2("echo -e '".implode("\n", $genes)."' | ".get_path("ngs-bits")."GenesToApproved");
	foreach($stdout as $line)
	{
		list($gene, $message) = explode("\t", $line);
		if (contains($message, "ERROR")) continue;
		if (in_array($gene, $genes_approved)) continue;
		$genes_approved[] = $gene;
	}
	
	//only with coordinates
	foreach($genes_approved as $gene)
	{
		list($stdout, $stderr, $exit_code) = exec2("echo '$gene' | ".get_path("ngs-bits")."GenesToBed -source ensembl -mode gene -fallback | ".get_path("ngs-bits")."BedMerge", false);
		if ($exit_code==0 && trim(implode("", $stdout))!="")
		{
			foreach($stdout as $line)
			{
				list($chr, $start, $end) = explode("\t", trim($line));
				if ($chr!=$chr_omim) continue;
				$start -= 20;
				$end += 20;
				print "$chr	$start	$end	{$mim_id}_[{$gene}_({$status})_{$disorders}]\n";
			}
		}
	}
}

fclose($handle);

?>