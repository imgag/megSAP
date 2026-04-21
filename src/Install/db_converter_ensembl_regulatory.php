<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//create ENSG>symbol map from Ensembl GFF file
//1       havana  ncRNA_gene      11121   24894   .       +       .       ID=gene:ENSG00000290825;Name=DDX11L16;biotype=lncRNA;description=DEAD/H-box helicase 11 like 16 (pseudogene) [Source:NCBI gene (formerly Entrezgene)%3BAcc:727856];gene_id=ENSG00000290825;logic_name=havana_homo_sapiens;version=2
//1       havana  pseudogene      12010   13670   .       +       .       ID=gene:ENSG00000223972;Name=DDX11L1;biotype=transcribed_unprocessed_pseudogene;description=DEAD/H-box helicase 11 like 1 (pseudogene) [Source:HGNC Symbol%3BAcc:HGNC:37102];gene_id=ENSG00000223972;logic_name=havana_homo_sapiens;version=6
$ensg2symbol = [];
$handle = gzopen2($argv[1], "r");
while (!feof($handle))
{
	//load
	$line = nl_trim(fgets($handle));
	if (!contains($line, "ID=gene:")) continue;
	
	$parts = explode("\t", $line);
	if (count($parts)<9) continue;
	
	$ensg = "";
	$symbol = "";
	foreach(explode(";", $parts[8]) as $entry)
	{
		if (starts_with($entry, "ID=gene:")) $ensg = substr($entry, 8);
		if (starts_with($entry, "Name=")) $symbol = substr($entry, 5);
	}
	if ($ensg!="" && $symbol!="")
	{
		$ensg2symbol[$ensg] = $symbol;
	}
}
fclose($handle);

//convert regulatory features
//1	Ensembl	promoter	10936	11436	.	.	.	ID=ENSR1_958;extended_start=9953;extended_end=11436;gene_id=ENSG00000290825;gene_name=DDX11L16;gene_biotype=lncRNA;color=#ff0000
//1	Ensembl	CTCF_binding_site	11222	11243	.	+	.	ID=ENSR1_53X2;color=#40e0d0
//1	Ensembl	enhancer	11437	11649	.	.	.	ID=ENSR1_88N;color=#faca00
$out = temp_file("_ctcf.bed");
$handle_out = fopen2($out, 'w');
$handle = gzopen2($argv[2], "r");
while (!feof($handle))
{
	//load
	$line = nl_trim(fgets($handle));
	$parts = explode("\t", $line);
	if (count($parts)<9) continue;
	
	$chr = "chr".trim($parts[0]);
	$type = trim($parts[2]);
	$start = trim($parts[3])-1;
	$end = trim($parts[4]);

	$id = "";
	$genes = [];
	foreach(explode(";", $parts[8]) as $entry)
	{
		if (starts_with($entry, "ID=")) $id = substr($entry, 3);
		if (starts_with($entry, "gene_id="))
		{
			$ensgs = explode(",", substr($entry, 8));
			foreach($ensgs as $ensg)
			{
				if (isset($ensg2symbol[$ensg])) $genes[] = $ensg2symbol[$ensg];
			}
		}
	}
	
	print implode("\t", [$chr, $start, $end, $id."|".$type."|".implode("&", $genes)])."\n";
	if ($type=="CTCF_binding_site")
	{
		fputs($handle_out, implode("\t", [$chr, $start, $end, $id."|".$type."|".implode("&", $genes)])."\n");
	}
}
fclose($handle);
fclose($handle_out);

//convert TF_binding_site motifs
//1	Ensembl	TF_binding_site	10144	10157	8.51609132380628	-	.	ID=ENSM00000000001;binding_matrix_id=ENSPFM0290;transcription_factor=HNF4A,NR2F1,RXRA
//1	Ensembl	TF_binding_site	10228	10241	8.51609132380628	-	.	ID=ENSM00000000002;binding_matrix_id=ENSPFM0290;transcription_factor=HNF4A,NR2F1,RXRA
$out2 = temp_file("_tf.bed");
$handle_out2 = fopen2($out2, 'w');
$handle = gzopen2($argv[3], "r");
while (!feof($handle))
{
	//load
	$line = nl_trim(fgets($handle));
	$parts = explode("\t", $line);
	if (count($parts)<9) continue;
	
	$chr = "chr".trim($parts[0]);
	$type = trim($parts[2]);
	$start = trim($parts[3])-1;
	$end = trim($parts[4]);
	
	fputs($handle_out2, implode("\t", [$chr, $start, $end, $type])."\n");
}
fclose($handle);
fclose($handle_out2);

//CTCF_binding_site are in both files, remove them from the TF_binding_site list
list($stdout) = exec2("BedSubtract -in $out2 -in2 $out | BedMerge");
foreach($stdout as $line)
{
	$line = nl_trim($line);
	if ($line=="") continue;
	print $line."\t|TF_binding_site|\n";
}
?>