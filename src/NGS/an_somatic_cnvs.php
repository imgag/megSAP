<?php
/** 
	@page an_somatic_cnvs
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_somatic_cnvs", "Annotates additional somatic data to ClinCNV file (NCG7.1 / RNA counts).");
$parser->addInfile("cnv_in", "Input CNV file.", false);
$parser->addFlag("include_ncg", "Annotate column with info from NCG7.1 whether a gene is TSG or oncogene");
$parser->addFlag("include_cytoband", "Annotate column with cytoband names.");
$parser->addInfile("rna_counts", "Annotate transcript counts from RNA counts file", true);
$parser->addString("rna_id", "Processed sample ID of RNA to be annotated.", true);
$parser->addString("rna_ref_tissue", "RNA tissue type for annotation of RNA reference transciptome data.", true);
$parser->addOutfile("out", "Output file name", false);
extract($parser->parse($argv));

if(!$include_ncg && !isset($rna_counts) && !$include_cytoband)
{
	trigger_error("At least one input file with annotation data must be given.", E_USER_ERROR);
}

$cnv_input = Matrix::fromTSV($cnv_in);

if($include_ncg)
{
	//remove old annotation
	$cnv_input->removeColByName("ncg_oncogene");
	$cnv_input->removeColByName("ncg_tsg");
	
	$handle = fopen2(get_path("data_folder") . "/dbs/NCG7.1/NCG7.1_oncogene.tsv", "r");
	$oncogenes =  array();
	$tsgs = array();
	while(!feof($handle))
	{
		$line = trim(fgets($handle));
		if(empty($line)) continue;
		if(starts_with($line, "entrez")) continue;
		
		list($entrezid,$gene, $cgc, $vogelstein, $oncogene, $tsg) = explode("\t", $line);
		if($oncogene == "1") $oncogenes[] = $gene;
		if($tsg == "1") $tsgs[] = $gene;
	}
	
	$oncogenes = approve_gene_names($oncogenes);
	$tsgs = approve_gene_names($tsgs);

	$col_tsg = array();
	$col_oncogene = array();
	$i_genes = $cnv_input->getColumnIndex("genes",false,false);
	
	if($i_genes === false) 
	{
		trigger_error("Cannot annotate file $cnv_in with NCG7.1 data because there is no column 'genes'.",E_USER_WARNING);
		exit(1);
	}
	
	//Annotate per gene
	for($row=0;$row<$cnv_input->rows();++$row)
	{
		$genes = explode(",",trim($cnv_input->get($row,$i_genes)));
		
		$col_oncogene[] = implode(",",array_intersect($genes, $oncogenes));
		$col_tsg[] = implode(",", array_intersect($genes, $tsgs));
	}
	$cnv_input->addCol($col_oncogene,"ncg_oncogene","oncogenes from NCG7.1 that overlap given CNV.");
	$cnv_input->addCol($col_tsg,"ncg_tsg","tumor suppressor genes from NCG7.1 that overlap given CNV.");
}

if($include_cytoband)
{
	//determine cytobands
	$cytobands = array();
	for($i=0;$i<$cnv_input->rows();++$i)
	{
		list($chr, $start, $end) = $cnv_input->getRow($i);
		$cytobands[] = implode(",", cytobands($chr, $start, $end));
	}
	$cnv_input->removeComment("cytoband", true);
	$cnv_input->removeColByName("cytoband");
	
	$cnv_input->addCol($cytobands, "cytoband", "Cytobands that are affected by CNV.");
}

if(isset($rna_counts))
{
	
	//Resubstitute zeroes by spaces (opposite happens in somatic_dna.php)
	$rna_ref_tissue = str_replace("0", " ", $rna_ref_tissue);
	
	if(!isset($rna_id) || !isset($rna_ref_tissue))
	{
		trigger_error("Both parameters \"-rna_ref_tissue\" and \"-rna_id\" have to be specified for the annotation of RNA data.", E_USER_ERROR);
	}
	
	//Make list of genes that occur in CNV file
	$i_genes = $cnv_input->getColumnIndex("genes", false, false);
	$genes_of_interest = array();
	for($row=0; $row<$cnv_input->rows(); ++$row)
	{
		$temp = explode( ",",trim($cnv_input->get($row, $i_genes)) );
		$genes_of_interest = array_merge($genes_of_interest, $temp);
	}
	$genes_of_interest = approve_gene_names(array_unique($genes_of_interest) );
	
	
	//create map of approved gene symbols in RNA counts file
	$approved_gene_map = array();
	$temp_file = temp_file(".tsv", "approved_gene_symbols");
	exec2("cut -f6 $rna_counts | sort | uniq | " . get_path("ngs-bits") . "/GenesToApproved | grep REPLACED > $temp_file", false);
	$handle = fopen2($temp_file, "r");
	while(!feof($handle))
	{
		$line = trim(fgets($handle));
		if(empty($line)) continue;
		
		list($new_symbol, $raw_string) = explode("\t", $line);
		
		//old gene symbol is contained as text in a sentence of the form "REPLACED: SYMBOL is a ..."
		list($old_symbol) = explode(" ", str_replace("REPLACED: ", "", $raw_string));
		$approved_gene_map[$old_symbol] = $new_symbol;
	}
	fclose($handle);
	
	//Create result array of genes and tpm that occur in RNA_counts and CNV file
	$results = array();
	$handle = fopen2($rna_counts, "r");
	while(!feof($handle))
	{
		$line = fgets($handle);
		
		if(starts_with($line,"#")) continue;
		if(empty($line)) continue;
		
		//order of line: gene_id, raw count, cpm, fpkm, tpm, gene symbol
		list(,,,,$tpm, $rna_gene) = explode("\t", $line);
		$rna_gene = strtoupper($rna_gene);
		
		//replace outdated gene symbols
		if(array_key_exists($rna_gene, $approved_gene_map)) 
		{
			$rna_gene = $approved_gene_map[$rna_gene];
		}
		
		if(in_array($rna_gene, $genes_of_interest))
		{
			$results[$rna_gene] = $tpm;
		}
	}
	
	//annotate RNA reference counts
	$ref_file = get_path("data_folder") . "/dbs/gene_expression/rna_tissue_consensus_v23.tsv";
	$ref_results = array();
	$handle = fopen2($ref_file,"r");
	$entry_count = 0;
	while(!feof($handle))
	{
		$line = fgets($handle);
		if(starts_with($line, "Gene\t")) continue;
		if(empty($line)) continue;

		list(,$ref_gene,$tissue,$tpm) = explode("\t", trim($line));
		
		if($tissue != trim($rna_ref_tissue)) continue;
		
		++$entry_count;
		
		if(in_array($ref_gene, $genes_of_interest)) $ref_results[$ref_gene] = $tpm;
	}
	fclose($handle);
	
	if($entry_count == 0)
	{
		trigger_error("Could not find any reference value for tissue type \"{$rna_ref_tissue}\".", E_USER_ERROR);
	}
	
	//Parse results to out file
	$new_col = array();
	
	//column that contains data from RNA reference
	$new_col_ref = array();
	
	for($row=0; $row<$cnv_input->rows(); ++$row)
	{
		$dna_genes = approve_gene_names(explode(",", trim($cnv_input->get($row, $i_genes))));
		$new_entry = array();
		$new_entry_ref =array();
		foreach($dna_genes as $dna_gene)
		{
			if(array_key_exists($dna_gene,$results))
			{
				$new_entry[] = "{$dna_gene}=". sprintf( "%.2f", round($results[$dna_gene], 2) );
			}
			else
			{
				$new_entry[] = "{$dna_gene}=.";
			}
			
			
			if(array_key_exists($dna_gene,$ref_results))
			{
				$new_entry_ref[] = "{$dna_gene}=". sprintf( "%.2f",round($ref_results[$dna_gene], 2) );
			}
			else
			{
				$new_entry_ref[] = "{$dna_gene}=.";
			}
		}
		$new_col[] = implode(",", $new_entry);
		$new_col_ref[] = implode(",", $new_entry_ref);
	}
	

	//remove outdated annotation columns/comments
	$cnv_input->removeComment("RNA_PROCESSED_SAMPLE_ID", true);
	$cnv_input->removeColByName("rna_tpm");
	$cnv_input->removeColByName("rna_ref_tpm");

	//Add annotation data
	$cnv_input->removeComment("RNA_REF_TPM_TISSUE=", true);
	$cnv_input->addComment("#RNA_REF_TPM_TISSUE={$rna_ref_tissue}");
	
	$cnv_input->removeColByName("{$rna_id}_rna_tpm");
	$cnv_input->addCol($new_col, "{$rna_id}_rna_tpm", "RNA TPM counts per gene from file {$rna_counts}.");
	
	$cnv_input->removeColByName("rna_ref_tpm", true);
	$cnv_input->addCol($new_col_ref, "rna_ref_tpm", "RNA reference data in transcripts per million for tissue {$rna_ref_tissue} from proteinatlas.org.");
}

$cnv_input->toTSV($out);

?>