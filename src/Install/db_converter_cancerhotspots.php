<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("db_converter_cancerhotspots", "Converts the CADD flat files (tsv.gz) to VCF(.gz).");
// optional
$parser->addInfile("in",  "Input file in tsv-format, copied from cancerhotspots.org xls-file", false, false);
$parser->addInfile("maf", "Input with cancerhotspots.org variants in maf-format.", false, false);
$parser->addOutfile("out",  "Output file in vcf.gz-format with CADD scores. ('-' for STDOUT)", true, "-");
extract($parser->parse($argv));


//get column indices of MAF file
$maf_handle = gzopen2($maf, "r");
while(!feof($maf_handle))
{
	$line = trim(fgets($maf_handle));
	if( starts_with($line, "Hugo_Symbol") )
	{
		$parts = explode("\t", $line);
		for($i=0; $i<count($parts); ++$i)
		{
			if($parts[$i] == "Hugo_Symbol") $i_maf_gene = $i;
			elseif($parts[$i] == "Entrez_Gene_Id") $i_maf_gene_id = $i;
			elseif($parts[$i] == "Transcript_ID") $i_maf_transcript = $i;
		}
		break;
	}
}
fclose($maf_handle);

// open input file
$handle = fopen2($in, "r");
// open output file
if ($out == "-") $out = "php://stdout";
$out_fp = fopen2($out, "w");

//write output header
fwrite($out_fp, "#GENE\tENTREZ_GENE_ID\tENSEMBL_TRANSCRIPT_ID\tAA_POS\tAA_REF\tAA_ALT\tCOUNT_TOTAL\tCOUNT_ALT\n");

while(!feof($handle))
{
	$line = trim(fgets($handle));
	
	//get column indices of input file
	if(starts_with($line, "Hugo_Symbol") )
	{
		$parts = explode("\t", $line);
		for($i=0; $i<count($parts); ++$i)
		{
			if($parts[$i] == "Hugo_Symbol") $i_gene = $i;
			elseif($parts[$i] == "Amino_Acid_Position") $i_aa_pos = $i;
			elseif($parts[$i] == "Mutation_Count") $i_mut_count = $i;
			elseif($parts[$i] == "Reference_Amino_Acid") $i_ref_aa = $i;
			elseif($parts[$i] == "Variant_Amino_Acid") $i_alt_aa = $i;
		}
		continue;
	}
	
	if(empty($line)) continue;
	
	$parts = explode("\t", $line);
	
	//get data from in to be included in output file
	$gene = $parts[$i_gene];
	$aa_pos = $parts[$i_aa_pos];
	$mut_count = $parts[$i_mut_count];
	$tmp_ref_aa = $parts[$i_ref_aa];
	$tmp_alt_aa = $parts[$i_alt_aa];
	
	//Parse reference and variant amino acid change
	list($ref_aa, $ref_count) = explode(":",$tmp_ref_aa);
	list($alt_aa, $alt_count) = explode(":",$tmp_alt_aa);
	
	if($ref_count != $mut_count)
	{
		trigger_error("invalid input data: $line \n", E_USER_ERROR);
	}


	//resolve transcript IDs from MAF file
	$maf_gene_id = "";
	$maf_transcript_id = "";
	
	$handle_maf = gzopen2($maf, "r");
	while(!feof($handle_maf))
	{
		$line_maf = trim(fgets($handle_maf));
		if(starts_with($line_maf, "#")) continue;
		if(starts_with($line_maf, "Hugo_Symbol")) continue;
		
		$parts = explode("\t", trim($line_maf));
		
		if($parts[$i_maf_gene] == $gene)
		{
			$maf_gene_id = $parts[$i_maf_gene_id];
			$maf_transcript_id = $parts[$i_maf_transcript];
			break;
		}
	}
	fclose($handle_maf);
	
	//write output line
	fwrite($out_fp, "{$gene}\t{$maf_gene_id}\t{$maf_transcript_id}\t{$aa_pos}\t{$ref_aa}\t{$alt_aa}\t{$mut_count}\t{$alt_count}\n");
}


fclose($handle);



?>