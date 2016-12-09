<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//load gene coordinates from UCSC
/*
#ID	chr	cds_start	cds_end	exon_starts	exon_ends	gene_name
uc010nxq.1	chr1	12189	13639	12189,12594,13402,	12227,12721,13639,	DDX11L1
uc001aal.1	chr1	69090	70008	69090,	70008,	OR4F5
uc021oeg.2	chr1	138529	139792	138529,139789,	139696,139792,	LOC729737
uc009vjk.2	chr1	324342	325605	324342,324438,	324345,325605,	LOC100133331
uc001aau.3	chr1	324342	325605	324342,324438,	324345,325605,	LOC100132062
uc021oeh.1	chr1	324514	325605	324514,324718,325382,	324686,325124,325605,	LOC100133331
*/
$gene2coord = array();
$file = file("../UCSC/kgXref_joined.txt");
foreach ($file as $line)
{
	if (trim($line)=="" || $line[0]=="#") continue;
	
	list(, $chr, $start, $end, , , $gene) = explode("\t", rtrim($line));
	
	//include only regular chromosomes
	if (contains($chr, "_")) continue;
	
	$hash = $chr."_".$gene;
	if (!isset($gene2coord[$hash]))
	{
		$gene2coord[$hash] = array($start, $end);
	}
	else
	{
		$gene2coord[$hash][0] = min($start, $gene2coord[$hash][0]);
		$gene2coord[$hash][1] = max($end, $gene2coord[$hash][1]);
	}
}

//parse "mim2genes.txt" to get approved symbols
/*
# Mim Number	Type	Gene IDs	Approved Gene Symbols
100050	predominantly phenotypes	-	-
100070	phenotype	100329167	-
100100	phenotype	-	-
100200	predominantly phenotypes	-	-
100300	phenotype	-	-
100500	moved/removed	-	-
100600	phenotype	-	-
100640	gene	216	ALDH1A1
100650	gene/phenotype	217	ALDH2
100660	gene	218	ALDH3A1
100670	gene	219	ALDH1B1
*/
$mim2genes = array();
$handle = fopen("mim2gene.txt", "r");
while(!feof($handle))
{
	$line = trim(fgets($handle));
	if ($line=="" || $line[0]=="#") continue;
	list($id, $type, $genes, $genes_approved) = explode("\t", $line);
	$genes_approved = trim($genes_approved);
	if ($genes_approved!="-")
	{
		$mim2genes[$id] = $genes_approved;
	}
}
fclose($handle);

//parse "genemap" for disorders
/*
1.1|5|13|13|1pter-p36.13|CTRCT8, CCV|P|Cataract, congenital, Volkmann type||115665|Fd|linked to Rh in Scottish family||Cataract 8, multiple types (2)| | ||
1.2|9|25|01|1p36.23|ENO1, PPH, MPB1|C|Enolase-1, alpha||172430|S, F, R, REa|||Enolase deficiency (1)| | |4(Eno1)|
1.3|12|22|87|1pter-p36|ERPL1, HLM2|C|Endogenous retroviral pol gene-like sequence 1 (oncogene HLM2)||131190|REa, F|||| | ||
1.4|4|14|11|1p36.11|HMGCL|P|3-hydroxy-3-methylglutaryl-Coenzyme A lyase||613898|REa, A|||HMG-CoA lyase deficiency, 246450 (3)| | |4(Hmgcl)|
1.5|4|29|14|1p36.33|AGRN, CMSPPD|P|Agrin||103320|REa|||Myasthenic syndrome, congenital, with pre- and postsynaptic|defects, 615120 (3) | |4(Agrn)|
1.6|3|15|92|1p36.33|GNB1|C|Guanine nucleotide-binding protein, beta polypeptide-1||139380|REa, A|||| | |4(Gnb1)|

For the file genemap, the fields are, in order :
1  - Numbering system, in the format  Chromosome.Map_Entry_Number
2  - Month entered
3  - Day     "
4  - Year    "
5  - Cytogenetic location
6  - Gene Symbol(s)
7  - Gene Status (see below for codes)
8  - Title
9  - Title, cont.
10 - MIM Number
11 - Method (see below for codes)
12 - Commentsq
13 - Comments, cont.
14 - Disorders (each disorder is followed by its MIM number, if
	different from that of the locus, and phenotype mapping method (see
	below).  Allelic disorders are separated by a semi-colon.
15 - Disorders, cont.
16 - Disorders, cont.
17 - Mouse correlate
18 - Reference
*/
$handle = fopen("genemap", "r");
print "#chr	start	end	omim\n";
while(!feof($handle))
{
	$line = fgets($handle);
	if(trim($line)=="")
	{
		continue;
	}
	
	list($chr, , , , , , $status, , , $mim_id, , , , $disorders1, $disorders2, $disorders3) = explode("|", $line);
	$mim_id = trim($mim_id);
	
	//only with approved symbol
	if (!isset($mim2genes[$mim_id])) continue;
	$genes = explode(",", $mim2genes[$mim_id]);
	
	//not limbo
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
	
	//disorders	
	$disorders = trim($disorders1." ".$disorders2." ".$disorders3);
	$disorders = preg_replace("/, \d{6,6}/i", "", $disorders); 
	$disorders = trim(strtr($disorders, array("(1)"=>"", "(2)"=>"", "(3)"=>"", "(4)"=>"", "["=>"", "]"=>"", "{"=>"", "}"=>"", "?"=>"", ";"=>"|")));
	$disorders = trim(strtr($disorders, array(" , "=>",", " ,"=>",", ", "=>",")));
	$disorders = trim(strtr($disorders, array(" | "=>"|", " |"=>"|", "| "=>"|")));
	$disorders = trim(strtr($disorders, array(" "=>"_")));
	if ($disorders=="") continue;
	$chr = "chr".substr($chr, 0, strpos($chr, "."));
	if ($chr=="chr23") $chr="chrX";
	if ($chr=="chr24") $chr="chrY";
	
	foreach($genes as $gene)
	{
		//only with coordinates
		if(!isset($gene2coord[$chr."_".$gene])) continue;
		list($start, $end) = $gene2coord[$chr."_".$gene];
		print "$chr	$start	$end	{$mim_id}_[{$gene}_({$status})_{$disorders}]\n";
	}
}

fclose($handle);

?>