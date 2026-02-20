<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("an_somatic_gsvar", "Annotates additional somatic data to GSvar file (NCG7.2 / RNA).");
$parser->addInfile("gsvar_in", "Input .gsvar-file with SNV data.", false);
$parser->addFlag("include_ncg", "Annotate column with info from NCG7.2 whether a gene is TSG or oncogene");
$parser->addString("rna_id", "ID of RNA sample if RNA data is annotated.", true);
$parser->addInfile("rna_counts", "Input file that contains RNA transcript counts.", true);
$parser->addInfile("rna_bam", "RNA-BAM file that is used to annotate and calculate variant depth and frequency in RNA sample.", true);
$parser->addString("rna_ref_tissue", "Annotate RNA reference data from The Human Protein Atlas in TPM (transcripts per million). Specify tissue, e.g. \"colon\". Zeroes will be replaced by blank space", true);
$parser->addString("rna_target", "BED file: limits the annotation of depth and frequency data to the given target. if 'NGSD' is given it tries to determine the target file from NGSD.", true, "");
$parser->addOutfile("out", "Output file name", false);
extract($parser->parse($argv));

$gsvar_input = Null;

if(isset($rna_bam) || isset($rna_counts) || isset($rna_id))
{
	if(!isset($rna_counts) || !isset($rna_id) || !isset($rna_bam))
	{
		trigger_error("For annotation of RNA data the parameters -rna_bam, -rna_counts and -rna_id must be specified jointly.", E_USER_ERROR);
	}
	
	$tmp_gsvar = $parser->tempFile(".Gsvar", "an_somatic_gsvar_");
	
	/************************
	 * RNA DEPTH AND RNA AF *
	 ************************/
	$parser->log("Annotation of RNA depth and af starting.");
	
	if ($rna_target == "NGSD") $rna_target = load_system($null, $rna_id)["target_file"];
	
	$args = [];
	$args[] = "-in {$gsvar_in}";
	$args[] = "-bam {$rna_bam}";
	$args[] = "-out {$tmp_gsvar}";
	$args[] = "-fragments";
	$args[] = "-depth";
	$args[] = "-name {$rna_id}";
	if ($rna_target != "") $args[] = "-target {$rna_target}";
	$parser->execApptainer("ngs-bits", "VariantAnnotateFrequency", implode(" ", $args), [$gsvar_in, $rna_bam, $rna_target], []);
	
	$gsvar_input = Matrix::fromTSV($tmp_gsvar);
	
	//Remove old annotations 
	$gsvar_input->removeColByName("{$rna_id}_rna_depth");
	$gsvar_input->removeColByName("{$rna_id}_rna_af");
	
	//Remove outdated comments/columns in header
	$gsvar_input->removeColByName("rna_depth");
	$gsvar_input->removeColByName("rna_af");
	$gsvar_input->removeComment("RNA_PROCESSED_SAMPLE_ID", true);
	$gsvar_input->removeComment("RNA_BAM_FILE", true);
	
	//rename headers and comments
	$gsvar_headers = $gsvar_input->getHeaders();
	
	$idx_freq = array_search("{$rna_id}_freq", $gsvar_headers);
	if ($idx_freq !== false) $gsvar_headers[$idx_freq] = "{$rna_id}_rna_af";
	
	$idx_depth = array_search("{$rna_id}_depth", $gsvar_headers);
	if ($idx_depth !== false) $gsvar_headers[$idx_depth] = "{$rna_id}_rna_depth";
	
	$gsvar_headers = $gsvar_input->setHeaders($gsvar_headers);
	
	$gsvar_comments = $gsvar_input->getComments();
	for($i=0; $i<count($gsvar_comments); $i++)
	{
		$comment = $gsvar_comments[$i];
		if (starts_with($comment, "#DESCRIPTION={$rna_id}_freq="))
		{
			$gsvar_comments[$i] = "#DESCRIPTION={$rna_id}_rna_af=Allele frequency in RNA BAM file ".realpath($rna_bam);
		}
		if (starts_with($comment, "#DESCRIPTION={$rna_id}_depth="))
		{
			$gsvar_comments[$i] = "#DESCRIPTION={$rna_id}_rna_depth=Depth in RNA BAM file ".realpath($rna_bam);
		}
	}
	$gsvar_input->setComments($gsvar_comments);
	
	$parser->log("Annotation of RNA depth and af finished.");
	
	/*********************
	 * TRANSCRIPT COUNTS *
	 *********************/
	$parser->log("Annotation of RNA transcription counts starting:");
	//Remove old annotations
	$gsvar_input->removeColByName("{$rna_id}_rna_tpm");
	$gsvar_input->removeColByName("rna_tpm");
	$gsvar_input->removeComment("RNA_TRANSCRIPT_COUNT_FILE", true);
	
	$handle = fopen2($rna_counts,"r");	
	
	$gene_to_tpm = array();
	$i_dna_gene =  $gsvar_input->getColumnIndex("gene");
	
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
		$rna_gene = "";
		if (contains($parts[$i_rna_gene], ","))
		{
			//use first gene in case there is more than one gene
			trigger_error("RNA count file contains multiple genes in gene_name column: ".$parts[$i_rna_gene], E_USER_NOTICE);
			$rna_gene = explode(",", $parts[$i_rna_gene])[0];
		}
		else
		{
			$rna_gene = $parts[$i_rna_gene];
		}
		if ($rna_gene == "") continue; //skip unnamed genes;
		$gene_to_tpm[$rna_gene] = number_format($parts[$i_rna_tpm], 2);
	}
	fclose($handle);
	
	
	$tpm_column = [];
	
	for($i=0; $i<$gsvar_input->rows(); ++$i)
	{	
		$genes = $gsvar_input->get($i,$i_dna_gene);
		
		if ($genes == "")
		{
			$tpm_column[] = ".";
			continue;
		}
		
		$genes = explode(",", $genes);
		
		$var_tpms = [];
		foreach($genes as $gene)
		{
			if (array_key_exists($gene, $gene_to_tpm)) $var_tpms[] = $gene_to_tpm[$gene];
			else $var_tpms[] = ".";
		}
		
		$tpm_column[] = implode("/", $var_tpms); 
	}
	
	$gsvar_input->addCol($tpm_column,"{$rna_id}_rna_tpm","Transcript count as annotated per gene from " . realpath($rna_counts));
	$parser->log("Annotation of RNA transcription counts finished.");
}

if ($gsvar_input ==  null) $gsvar_input = Matrix::fromTSV($gsvar_in);

if(isset($rna_ref_tissue))
{
	$parser->log("Annotation of reference tissue data starting:");
	
	//Resubstitute zeroes by spaces (opposite happens in somatic_tumor_normal.php/somatic_tumor_only.php)
	$rna_ref_tissue = str_replace("0", " ", $rna_ref_tissue);
	
	$ref_file = get_path("data_folder") . "/dbs/gene_expression/rna_tissue_consensus_v24.tsv";
	
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

	$gene_to_ref_tpm = [];
	
	while(!feof($handle))
	{
		$line = trim(fgets($handle));
		if(empty($line)) continue; 
		$parts = explode("\t",$line);
		
		if($parts[$i_ref_tissue_type] != $rna_ref_tissue) continue;
		
		$gene_to_ref_tpm[$parts[$i_ref_gene]] = $parts[$i_ref_tpm];
	}
	fclose($handle);
	
	if(count($gene_to_ref_tpm) == 0)
	{
		trigger_error("Could not find any entry for tissue type $rna_ref_tissue in {$ref_file}.", E_USER_ERROR);
	}
	
	
	$i_dna_gene =  $gsvar_input->getColumnIndex("gene");
	$ref_tpm_column = [];
	
	for($i=0; $i<$gsvar_input->rows(); ++$i)
	{	
		$genes = $gsvar_input->get($i,$i_dna_gene);
		
		if ($genes == "")
		{
			$ref_tpm_column[] = ".";
			continue;
		}
		
		$genes = explode(",", $genes);
		
		$var_tpms = [];
		foreach($genes as $gene)
		{
			if (array_key_exists($gene, $gene_to_ref_tpm)) $var_tpms[] = $gene_to_ref_tpm[$gene];
			else $var_tpms[] = ".";
		}
		
		$ref_tpm_column[] = implode("/", $var_tpms); 
	}
	
	$gsvar_input->removeComment("RNA_REF_TPM_TISSUE=", true);
	$gsvar_input->addComment("#RNA_REF_TPM_TISSUE={$rna_ref_tissue}");
	
	$gsvar_input->removeColByName("rna_ref_tpm");
	$gsvar_input->addCol($ref_tpm_column, "rna_ref_tpm", "RNA reference data in transcripts per million for tissue {$rna_ref_tissue} from proteinatlas.org.");
	$parser->log("Annotation of reference tissue data finished.");
}

if($include_ncg)
{
	$parser->log("Annotation of NCG starting:");
	$ncg_file = get_path("data_folder") . "/dbs/NCG7.2/NCG7.2_oncogene.tsv";
	annotate_gsvar_by_gene($gsvar_input, $ncg_file, "symbol", "NCG_oncogene", "ncg_oncogene", "1:gene is oncogene according NCG7.2, 0:No oncogene according NCG7.2, na: no information available about gene in NCG7.2. Order is the same as in column gene.", false);
	annotate_gsvar_by_gene($gsvar_input, $ncg_file, "symbol", "NCG_tsg", "ncg_tsg", "1:gene is TSG according NCG7.2, 0:No TSG according NCG7.2, na: no information available about gene in NCG7.2. Order is the same as in column gene.", false);	
	$parser->log("Annotation of NCG finished.");
}

$gsvar_input->toTSV($out);

?>
