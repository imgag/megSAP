<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//init
$ngsbits = get_path("ngs-bits");
$c_chr_mismatch = 0;
$c_no_approved_gene = 0;

//parse "genemap2.txt" for coordinates and disorders
/*
# Chromosome	Genomic Position Start	Genomic Position End	Cyto Location	Computed Cyto Location	Mim Number	Gene Symbols	Gene Name	Approved Symbol	Entrez Gene ID	Ensembl Gene ID	Comments	Phenotypes	Mouse Gene Symbol/ID
chr1	975197	982116	1p36.33	1p36.33	615921	PERM1, C1orf170	PPARGC1-and ESRR-induced regulator, muscle, 1	PERM1	84808	ENSG00000187642			
chr1	998961	1001051	1p36.31	1p36.33	608060	HES4	Hairy/enhancer of split, Drosophila, homolog of, 4	HES4	57801	ENSG00000188290			
*/
$handle = fopen2("genemap2.txt", "r");
print "#chr	start	end	omim\n";
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if($line=="" || $line[0]=="#") continue;
	
	$parts = explode("\t", $line);
	if (count($parts)<14)
	{
		print_r($parts);
		trigger_error("Line '$line' does not contain 14 tab-separated parts!", E_USER_ERROR);
	}
	
	//id
	$mim_id = trim($parts[5]);
	
	//disorders not empty
	$disorders = trim($parts[12]);
	$disorders = trim(strtr($disorders, array("(1)"=>"", "(2)"=>"", "(3)"=>"", "(4)"=>"", "["=>"", "]"=>"", "{"=>"", "}"=>"", "?"=>"", ";"=>"|")));
	$disorders = trim(strtr($disorders, array(" , "=>",", " ,"=>",", ", "=>",")));
	$disorders = trim(strtr($disorders, array(" | "=>"/", " |"=>"/", "| "=>"/")));
	$disorders = trim(strtr($disorders, array(" "=>"_")));
	if ($disorders=="") continue;
	
	//chromosome
	$chr_omim = trim($parts[0]);
	
	//extract genes
	$genes = array();
	foreach(explode(",", $parts[6].",".$parts[8]) as $gene)
	{
		$gene = trim($gene);
		if ($gene=="") continue;
		
		$genes[] = $gene;
	}
	
	//convert genes to approved symbols
	$genes_approved = array();
	list($stdout) = exec2("echo -e '".implode("\n", $genes)."' | {$ngsbits}GenesToApproved");
	foreach($stdout as $line)
	{
		list($gene, $message) = explode("\t", $line);
		if (contains($message, "ERROR")) continue;
		if (in_array($gene, $genes_approved)) continue;
		$genes_approved[] = $gene;
	}
	if (count($genes_approved)==0)
	{
		//print "##Warning: No approved genes for MIM {$mim_id}!\n";
		++$c_no_approved_gene;
		continue;
	}
	
	//only with coordinates
	foreach($genes_approved as $gene)
	{
		list($stdout, $stderr, $exit_code) = exec2("echo '$gene' | {$ngsbits}GenesToBed -source ensembl -mode gene -fallback | {$ngsbits}BedExtend -n 20 | {$ngsbits}BedMerge", false);
	
		if ($exit_code==0 && trim(implode("", $stdout))!="")
		{
			foreach($stdout as $line)
			{
				list($chr, $start, $end) = explode("\t", trim($line));
				if ($chr!=$chr_omim)
				{
					//print "##Warning: For MIM/gene '$mim_id/$gene' chromosome is '$chr', but OMIM chromosome is '$chr_omim'!\n";
					++$c_chr_mismatch;
					continue;
				}
				print "$chr	$start	$end	{$mim_id}_[GENE={$gene}_PHENOS={$disorders}]\n";
			}
		}
	}
}
fclose($handle);

print "##MIM entries without approved gene symbols: {$c_no_approved_gene}\n";
print "##MIM entries with chromosome mismatch: {$c_chr_mismatch}\n";

?>