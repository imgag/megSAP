<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("an_somatic_gsvar", "Annotates additional somatic data to GSvar file (NCG7.1 / RNA).");
$parser->addInfile("gsvar_in", "Input .gsvar-file with SNV data.", false);
$parser->addFlag("include_ncg", "Annotate column with info from NCG7.1 whether a gene is TSG or oncogene");
$parser->addString("rna_id", "ID of RNA sample if RNA data is annotated.", true);
$parser->addInfile("rna_counts", "Input file that contains RNA transcript counts.", true);
$parser->addInfile("rna_bam", "RNA-BAM file that is used to annotate and calculate variant depth and frequency in RNA sample.", true);
$parser->addString("rna_ref_tissue", "Annotate RNA reference data from The Human Protein Atlas in TPM (transcripts per million). Specify tissue, e.g. \"colon\". Zeroes will be replaced by blank space", true);
$parser->addOutfile("out", "Output file name", false);
extract($parser->parse($argv));

$gsvar_input = Matrix::fromTSV($gsvar_in);

if(isset($rna_bam) || isset($rna_counts) || isset($rna_id))
{
	if(!isset($rna_counts) || !isset($rna_id) || !isset($rna_bam))
	{
		trigger_error("For annotation of RNA data the parameters -rna_bam, -rna_counts and -rna_id must be specified jointly.", E_USER_ERROR);
	}
	
	/************************
	 * RNA DEPTH AND RNA AF *
	 ************************/
	//Remove old annotations 
	$gsvar_input->removeColByName("{$rna_id}_rna_tpm");
	$gsvar_input->removeColByName("{$rna_id}_rna_depth");
	$gsvar_input->removeColByName("{$rna_id}_rna_af");
	
	//Remove outdated comments/columns in header
	$gsvar_input->removeColByName("rna_tpm");
	$gsvar_input->removeColByName("rna_depth");
	$gsvar_input->removeColByName("rna_af");
	$gsvar_input->removeComment("RNA_PROCESSED_SAMPLE_ID", true);
	$gsvar_input->removeComment("RNA_TRANSCRIPT_COUNT_FILE", true);
	$gsvar_input->removeComment("RNA_BAM_FILE", true);
	 
	 
	$rna_afs = array();
	$rna_depths = array();
	
	$i_variant_type = $gsvar_input->getColumnIndex("variant_type");
	
	for($i=0; $i<$gsvar_input->rows(); ++$i)
	{
		//Skip intronic variants (in these cases we only see DNA artefacts)
		$variant_type =  $gsvar_input->get($i, $i_variant_type);
		if($variant_type == "intron" || $variant_type == "intergenic")
		{
			$rna_depths[] = ".";
			$rna_afs[] = ".";
			continue;
		}
		
		list($chr,$start,$end,$ref,$obs) = $gsvar_input->getRow($i);

		$counts = allele_count($rna_bam,$chr,$start);
		
		$depth = 0; //we use depth that only counts bases from mpileup and e.g. no reference skips
		foreach(array_values($counts) as $val) $depth += $val; 
		
		$obs = str_replace("-", "*", $obs); //deletion
		$obs = str_replace("+", "*", $obs); //insertion
		if(strlen($obs) != 1) $obs = "*"; //insertion 
		
		$rna_depths[] = $depth;
		$rna_afs[] = ($depth != 0) ? number_format($counts[$obs] / $depth,  2) : number_format(0, 2);
		
		if(!array_key_exists($obs,$counts)) echo $obs ."\n";
	}
	
	$gsvar_input->addCol($rna_depths, "{$rna_id}_rna_depth", "Depth in RNA BAM file " . realpath($rna_bam));
	$gsvar_input->addCol($rna_afs, "{$rna_id}_rna_af", "Allele frequency in RNA BAM file " . realpath($rna_bam));

	
	/*********************
	 * TRANSCRIPT COUNTS *
	 *********************/
	$handle = fopen2($rna_counts,"r");	
	
	$genes_of_interest = array();
	$i_dna_gene =  $gsvar_input->getColumnIndex("gene");
	
	//Make list of genes that appear in GSvar file, as associative array
	for($i=0; $i<$gsvar_input->rows(); ++$i)
	{
		list($chr,$start,$end,$ref,$obs) = $gsvar_input->getRow($i);
		
		$genes = explode(",", $gsvar_input->get($i,$i_dna_gene));
		
		$genes_of_interest["{$chr}_{$start}_{$end}_{$ref}_{$obs}"] = $genes;
	}
	
	$i_rna_gene = -1;
	$i_rna_tpm = -1;
	
	$results  = array();
	
	$db_is_enabled = db_is_enabled("NGSD");
	
	while(!feof($handle))
	{
		$line = trim(fgets($handle));
		
		//get indices of RNA gene and RNA TPM
		if(starts_with($line,"#gene_id"))
		{
			$parts = explode("\t", $line);
			for($i=0; $i<count($parts); ++$i)
			{
				if($parts[$i] ==  "gene_name") $i_rna_gene = $i;
				if($parts[$i] == "tpm") $i_rna_tpm = $i;
			}
		}
		
		//skip comments
		if(starts_with($line,"#")) continue;
		if(empty($line)) continue;
		
		
		$parts = explode("\t", $line);
		
		if($db_is_enabled) list($rna_gene) = (explode(",", $parts[$i_rna_gene])); //use first gene in case there is more than one gene
		else list($rna_gene) = explode(",", $parts[$i_rna_gene]);
		
		foreach($genes_of_interest as $key => $dna_genes)
		{
			if(!array_key_exists($key, $results)) $results[$key] = "."; //Make a dummy entry "." in case we find nothing
			
			
			foreach($dna_genes as $dna_gene)
			{
				if($dna_gene == $rna_gene)
				{
					$results[$key] = number_format($parts[$i_rna_tpm], 2);
				}
			}
		}
	}
	fclose($handle);
	
	$gsvar_input->addCol(array_values($results),"{$rna_id}_rna_tpm","Transcript count as annotated per gene from " . realpath($rna_counts));
}

if(isset($rna_ref_tissue))
{
	//Resubstitute zeroes by spaces (opposite happens in somatic_dna.php)
	$rna_ref_tissue = str_replace("0", " ", $rna_ref_tissue);
	
	$genes_of_interest = array();
	$i_dna_gene =  $gsvar_input->getColumnIndex("gene");
	
	//Make list of genes that appear in GSvar file, as associative array
	for($i=0; $i<$gsvar_input->rows(); ++$i)
	{
		list($chr,$start,$end,$ref,$obs) = $gsvar_input->getRow($i);
		$genes_of_interest["{$chr}_{$start}_{$end}_{$ref}_{$obs}"] = $gsvar_input->get($i,$i_dna_gene);
	}
	
	
	$ref_file = get_path("data_folder") . "/dbs/gene_expression/rna_tissue_consensus_v23.tsv";
	
	if(!file_exists($ref_file))
	{
		trigger_error("Could not find $ref_file neccessary for annotation of reference RNA expression.", E_USER_ERROR);
	}
	$handle = fopen2($ref_file, "r");
	
	$i_ref_gene = -1;
	$i_ref_tpm = -1;
	$i_ref_tissue_type = -1;
	
	//first line of $ref_file contains header
	$header_parts = explode("\t", fgets($handle));
	for($i=0; $i<count($header_parts); ++$i)
	{
		$header = trim($header_parts[$i]);
		if($header == "Gene name") $i_ref_gene = $i;
		if($header == "nTPM") $i_ref_tpm = $i;
		if($header == "Tissue") $i_ref_tissue_type = $i;
	}

	//ref_results contains contents of column with ref values to be annotated
	$ref_results = array_fill_keys(array_keys($genes_of_interest), ".");
	
	$ref_entry_count = 0;
	while(!feof($handle))
	{
		$line = trim(fgets($handle));
		if(empty($line)) continue; 
		$parts = explode("\t",$line);
		
		if($parts[$i_ref_tissue_type] != $rna_ref_tissue) continue; 
		
		foreach($genes_of_interest as $key => $dna_genes)
		{
			$expression_data = array();
			
			$tmp_dna_genes = explode(",", $dna_genes);
			foreach($tmp_dna_genes as $dna_gene)
			{
				if($dna_gene == $parts[$i_ref_gene])
				{
					$expression_data[] = $parts[$i_ref_tpm];
				}
			}
			if(!empty($expression_data))$ref_results[$key] = implode(",", $expression_data);
		}
		
		++$ref_entry_count;
	}
	
	if($ref_entry_count == 0)
	{
		trigger_error("Could not find any entry for tissue type $rna_ref_tissue in {$ref_file}.", E_USER_ERROR);
	}
	
	$gsvar_input->removeComment("RNA_REF_TPM_TISSUE=", true);
	$gsvar_input->addComment("#RNA_REF_TPM_TISSUE={$rna_ref_tissue}");
	
	$gsvar_input->removeColByName("rna_ref_tpm");
	$gsvar_input->addCol(array_values($ref_results), "rna_ref_tpm", "RNA reference data in transcripts per million for tissue {$rna_ref_tissue} from proteinatlas.org.");
}

if($include_ncg)
{
	$ncg_file = get_path("data_folder") . "/dbs/NCG7.1/NCG7.1_oncogene.tsv";
	annotate_gsvar_by_gene($gsvar_input, $ncg_file, "symbol", "NCG_oncogene", "ncg_oncogene", "1:gene is oncogene according NCG7.1, 0:No oncogene according NCG7.1, na: no information available about gene in NCG7.1. Order is the same as in column gene.", false);
	annotate_gsvar_by_gene($gsvar_input, $ncg_file, "symbol", "NCG_tsg", "ncg_tsg", "1:gene is TSG according NCG7.1, 0:No TSG according NCG7.1, na: no information available about gene in NCG7.1. Order is the same as in column gene.", false);	
}


$gsvar_input->toTSV($out);

?>
