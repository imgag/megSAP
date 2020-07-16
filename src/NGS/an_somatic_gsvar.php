<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("an_somatic_gsvar", "Annotates additional somatic data to GSvar file (CGI / NCG6.0).");
$parser->addInfile("gsvar_in", "Input .gsvar-file with SNV data.", false);
$parser->addInfile("cgi_snv_in", "Input CGI data with SNV annotations",true);
$parser->addFlag("cgi_empty", "Add empty CGI columns for convenience e.g. if there is no CGI output for variants of sample.");
$parser->addString("cgi_empty_cancer_type", "Cancer type for empty CGI columns.", true, "CANCER");
$parser->addFlag("include_ncg", "Annotate column with info from NCG6.0 whether a gene is TSG or oncogene");
$parser->addString("rna_id", "ID of RNA sample if RNA data is annotated.", true);
$parser->addInfile("rna_counts", "Input file that contains RNA transcript counts.", true);
$parser->addInfile("rna_bam", "RNA-BAM file that is used to annotate and calculate variant depth and frequency in RNA sample.", true);
$parser->addString("rna_ref_tissue", "Annotate RNA reference data from The Human Protein Atlas in TPM (transcripts per million). Specify tissue, e.g. \"colon\". Zeroes will be replaced by blank space", true);
$parser->addOutfile("out", "Output file name", false);
extract($parser->parse($argv));

//transforms CGI mutation file to VCF and sorts chrom. positions
//CHROM	POS	REF	ALT	CGI
function cgi_variant_tsv_to_sorted_vcf($cgi_input)
{
	global $parser;
	
	$cgi_file = Matrix::fromTSV($cgi_input);
	$output_lines = array("#CHROM\tSTART\tEND\tREF\tALT\tCGI\n");
	
	//create array containing CGI information
	//id,driver,gene_role,transcript,gene,consequence
	$i_id = $cgi_file->getColumnIndex("sample",false,false); //contains ID for each passed variant;	
	$i_drv = $cgi_file->getColumnIndex("driver_statement");
	$i_role = $cgi_file->getColumnIndex("gene_role");
	$i_trans = $cgi_file->getColumnIndex("transcript");
	$i_gene = $cgi_file->getColumnIndex("gene");
	$i_cons = $cgi_file->getColumnIndex("consequence");
	$cgi_statements = array();
	for($i=0;$i<$cgi_file->rows();++$i)
	{
		list($chr,$start,$end,$ref,$alt) =  explode( "_", $cgi_file->get($i, $i_id) ); //genomic coordinates extracted from id in the form "chr_start_end_ref_alt"
		
		$cgi_data = [
			str_replace(",",";",trim($cgi_file->get($i,$i_id))),
			str_replace(",",";",trim($cgi_file->get($i,$i_drv))),
			(empty(trim($cgi_file->get($i,$i_role))) ? "." : str_replace(",",";",trim($cgi_file->get($i,$i_role))) ),
			str_replace(",",";",trim($cgi_file->get($i,$i_trans))),
			str_replace(",",";",trim($cgi_file->get($i,$i_gene))),
			str_replace(",",";",trim($cgi_file->get($i,$i_cons)))
		];
		
		$output_lines[] = "{$chr}\t{$start}\t{$end}\t{$ref}\t{$alt}\t". implode(",",$cgi_data) ."\n";
	}
	
	$output_lines = array_unique($output_lines);
	
	//save to temp file and use ngs-bits for sorting
	$temp_o_file = $parser->tempFile(".tsv");
	file_put_contents($temp_o_file,$output_lines);
	$output = $parser->tempFile(".tsv");
	
	//sort file
	exec2("( head -n1 {$temp_o_file}; tail -n+2 {$temp_o_file} | sort -k1.4,1 -k2,2 -n )   > {$output} ");
	return $output;
}

/********
 * MAIN *
 ********/
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
	$handle = fopen($rna_counts,"r");	
	
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
	
	
	$ref_file = get_path("data_folder") . "/dbs/gene_expression/rna_tissue_hpa.tsv";
	
	$handle = fopen($ref_file, "r");
	
	$i_ref_gene = -1;
	$i_ref_tpm = -1;
	$i_ref_tissue_type = -1;
	
	//first line of $ref_file contains header
	$header_parts = explode("\t", fgets($handle));
	for($i=0; $i<count($header_parts); ++$i)
	{
		if($header_parts[$i] == "Gene name") $i_ref_gene = $i;
		if($header_parts[$i] == "TPM") $i_ref_tpm = $i;
		if($header_parts[$i] == "Tissue") $i_ref_tissue_type = $i;
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

if(isset($cgi_snv_in))
{
	/******************************
	 * REMOVE OLD CGI ANNOTATIONS *
	 ******************************/ 
	//old annotations created by former script
	$old_annotations = ["CGI_drug_assoc","CGI_evid_level","CGI_transcript","CGI_id","CGI_driver_statement","CGI_gene_role","CGI_gene","CGI_consequence"];
	foreach($old_annotations as $annotation_name) $gsvar_input->removeColByName($annotation_name);

	/***************************
	 * REMOVE OLD CGI COMMENTS *
	 ***************************/
	$old_comments = ["CGI_ICD10_CODE","CGI_ICD10_TEXT","CGI_ICD10_DIAGNOSES","GENES_FOR_REIMBURSEMENT","CGI_CANCER_TYPE","CGI_id","CGI_driver_statement","CGI_gene_role","CGI_transcript","CGI_gene","CGI_consequence"];
	foreach($old_comments as $old_comment) $gsvar_input->removeComment($old_comment,true);

	
	//VCF file with CGI mutation analysis
	$cgi_snvs = file(cgi_variant_tsv_to_sorted_vcf($cgi_snv_in),FILE_IGNORE_NEW_LINES);
	
	//new columns that will be filled using CGI mutation file
	$col_cgi_id = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_driver_statement = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_gene_role = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_transcript = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_gene = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_consequence = array_fill(0,$gsvar_input->rows(),"");
	
	
	//annotate GSvar file with CGI data, use positions of SNV for annotation
	for($i=0;$i<$gsvar_input->rows();$i++)
	{
		$gsvar_chr = $gsvar_input->get($i,0);
		$gsvar_start = $gsvar_input->get($i,1);
		$gsvar_end = $gsvar_input->get($i,2);
		$gsvar_ref = $gsvar_input->get($i,3);
		$gsvar_obs = $gsvar_input->get($i,4);
		
		foreach($cgi_snvs as $cgi_line)
		{
			
			if(starts_with($cgi_line,"#")) continue;
			list($cgi_chr,$cgi_start,$cgi_end,$cgi_ref,$cgi_obs,$cgi_statements) = explode("\t",$cgi_line);

			if((is_numeric($cgi_chr) && is_numeric($gsvar_chr)) && ($cgi_chr > $gsvar_chr)) break;
			if($cgi_chr != $gsvar_chr) continue;
			if($gsvar_start != $cgi_start) continue;
			if($gsvar_end != $cgi_end) continue;
			if($gsvar_ref != $cgi_ref) continue;
			if($gsvar_obs != $cgi_obs) continue;

			$cgi_statements = explode(",",$cgi_statements);
			
			$col_cgi_id[$i] = $cgi_statements[0];
			$col_cgi_driver_statement[$i] = $cgi_statements[1];
			$col_cgi_gene_role[$i] = $cgi_statements[2];
			$col_cgi_transcript[$i] = $cgi_statements[3];
			$col_cgi_gene[$i] = $cgi_statements[4];
			$col_cgi_consequence[$i] = $cgi_statements[5];
			
			break;
		}
	}
	//get cancer type from CGI input file
	$cgi_file = Matrix::fromTSV($cgi_snv_in);
	$cancer_type_cgi = $cgi_file->get(0,$cgi_file->getColumnIndex("cancer"));

	//add CGI cancer type information to GSVar
	if(!empty($cancer_type_cgi))	$gsvar_input->addComment("#CGI_CANCER_TYPE=$cancer_type_cgi");

	//save results to output GSvar file
	$gsvar_input->addCol($col_cgi_id,"CGI_id","Identifier for CGI statements");
	$gsvar_input->addCol($col_cgi_driver_statement,"CGI_driver_statement","CancerGenomeInterpreter.org oncogenic classification");
	$gsvar_input->addCol($col_cgi_gene_role,"CGI_gene_role","CancerGenomeInterpreter.org gene role. LoF: Loss of Function, Act: Activating");
	$gsvar_input->addCol($col_cgi_transcript,"CGI_transcript","CancerGenomeInterpreter.org CGI Ensembl transcript ID");
	$gsvar_input->addCol($col_cgi_gene,"CGI_gene","Gene symbol returned by CancerGenomeInterpreter.org");
	$gsvar_input->addCol($col_cgi_consequence,"CGI_consequence","Consequence of the mutation assessed by CancerGenomeInterpreter.org");
}

//Add empty CGI columns for convenience
if($cgi_empty)
{
	/******************************
	 * REMOVE OLD CGI ANNOTATIONS *
	 ******************************/ 
	//old annotations created by former script
	$old_annotations = ["CGI_drug_assoc","CGI_evid_level","CGI_transcript","CGI_id","CGI_driver_statement","CGI_gene_role","CGI_gene","CGI_consequence"];
	foreach($old_annotations as $annotation_name) $gsvar_input->removeColByName($annotation_name);

	/***************************
	 * REMOVE OLD CGI COMMENTS *
	 ***************************/
	$old_comments = ["CGI_ICD10_CODE","CGI_ICD10_TEXT","CGI_ICD10_DIAGNOSES","GENES_FOR_REIMBURSEMENT","CGI_CANCER_TYPE","CGI_id","CGI_driver_statement","CGI_gene_role","CGI_transcript","CGI_gene","CGI_consequence"];
	foreach($old_comments as $old_comment) $gsvar_input->removeComment($old_comment,true);
	
	
	$empty_col = array_fill(0, $gsvar_input->rows(), "");
	
	$gsvar_input->addCol($empty_col,"CGI_id","Identifier for CGI statements");
	$gsvar_input->addCol($empty_col,"CGI_driver_statement","CancerGenomeInterpreter.org oncogenic classification");
	$gsvar_input->addCol($empty_col,"CGI_gene_role","CancerGenomeInterpreter.org gene role. LoF: Loss of Function, Act: Activating");
	$gsvar_input->addCol($empty_col,"CGI_transcript","CancerGenomeInterpreter.org CGI Ensembl transcript ID");
	$gsvar_input->addCol($empty_col,"CGI_gene","Gene symbol returned by CancerGenomeInterpreter.org");
	$gsvar_input->addCol($empty_col,"CGI_consequence","Consequence of the mutation assessed by CancerGenomeInterpreter.org");
	
	$gsvar_input->addComment("#CGI_CANCER_TYPE=$cgi_empty_cancer_type");
}

if($include_ncg)
{
	//Remove potential old NCG annotation
	$gsvar_input->removeColByName("ncg_oncogene");
	$gsvar_input->removeColByName("ncg_tsg");
	
	$col_oncogene  = array();
	$col_tsg  = array();
	
	for($row=0;$row<$gsvar_input->rows();++$row)
	{
		//Get gene names in GSvar file
		$genes = explode(",",$gsvar_input->get($row,$gsvar_input->getColumnIndex("gene")));
		
		$ncg_oncogene = "";
		$ncg_tsg = "";
		
		//Annotate NCG information per gene
		foreach($genes as $gene)
		{
			$statements = ncg_gene_statements($gene);
			
			$ncg_oncogene .= $statements["is_oncogene"] .",";
			$ncg_tsg .= $statements["is_tsg"] . ",";
		}
		$ncg_oncogene = substr($ncg_oncogene,0,-1);
		$ncg_tsg = substr($ncg_tsg,0,-1);
		
		$col_oncogene[] = $ncg_oncogene;
		$col_tsg[] = $ncg_tsg;
	}
	
	$gsvar_input->addCol($col_oncogene,"ncg_oncogene","1:gene is oncogene according NCG6.0, 0:No oncogene according NCG6.0, na: no information available about gene in NCG6.0. Order is the same as in column gene.");
	$gsvar_input->addCol($col_tsg,"ncg_tsg","1:gene is TSG according NCG6.0, 0:No TSG according NCG6.0, na: no information available about gene in NCG6.0. Order is the same as in column gene.");
}


$gsvar_input->toTSV($out);

?>
