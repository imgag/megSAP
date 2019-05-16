<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("an_somatic_gsvar", "Annotates additional somatic data to GSvar file (CGI / NCG6.0).");
$parser->addInfile("gsvar_in", "Input .gsvar-file with SNV data.", false);
$parser->addInfile("cgi_snv_in", "Input CGI data with SNV annotations",true);
$parser->addFlag("include_ncg", "Annotate column with info from NCG6.0 whether a gene is TSG or oncogene");
$parser->addOutfile("out", "Output file name", false);
extract($parser->parse($argv));

//transforms CGI mutation file to VCF and sorts chrom. positions
//CHROM	POS	REF	ALT	CGI
function cgi_variant_tsv_to_sorted_vcf($cgi_input)
{
	$cgi_file = Matrix::fromTSV($cgi_input);
	$output_lines = array("#CHROM\tPOS\tREF\tALT\tCGI\n");
	
	//Chromosomal positions
	$i_pos = $cgi_file->getColumnIndex("pos",false,false);
	$i_chr = $cgi_file->getColumnIndex("chr",false,false);
	$i_ref = $cgi_file->getColumnIndex("ref",false,false);
	$i_alt = $cgi_file->getColumnIndex("alt",false,false);
	
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
		$temp = [
			str_replace(",",";",trim($cgi_file->get($i,$i_id))),
			str_replace(",",";",trim($cgi_file->get($i,$i_drv))),
			(empty(trim($cgi_file->get($i,$i_role))) ? "." : str_replace(",",";",trim($cgi_file->get($i,$i_role))) ),
			str_replace(",",";",trim($cgi_file->get($i,$i_trans))),
			str_replace(",",";",trim($cgi_file->get($i,$i_gene))),
			str_replace(",",";",trim($cgi_file->get($i,$i_cons)))
		];
		$cgi_statements[] = implode(",",$temp);
	}
	
	if($i_pos !== false && $i_chr !== false && $i_ref !== false && $i_alt !== false) //chr. positions included in file
	{
		for($i=0;$i<$cgi_file->rows();++$i)
		{
			$chr = "chr" . $cgi_file->get($i,$i_chr);
			$pos = $cgi_file->get($i,$i_pos);
			$ref = $cgi_file->get($i,$i_ref);
			$alt = $cgi_file->get($i,$i_alt);
			$output_lines[] = "{$chr}\t{$pos}\t{$ref}\t{$alt}\t". $cgi_statements[$i] ."\n";
		}
	}
	else //chrom. pos. can be found in input column if pos. do not have own columns
	{
		$i_input = $cgi_file->getColumnIndex("input",false,false);
		if($i_input === false)
		{
			trigger_error("Could not determine chromosomal positions in CGI file {$cgi_snv_in}. Aborting.",E_USER_ERROR);
		}
		$input = $cgi_file->getCol($i_input);
		for($i=0;$i<$cgi_file->rows();++$i)
		{
			list($chr,$pos,$ref,$alt) = explode("|",$input[$i]);
			$output_lines[] = "chr{$chr}\t{$pos}\t{$ref}\t{$alt}\t". $cgi_statements[$i] ."\n";
		}
	}
	
	//save to temp file and use ngs-bits for sorting
	$temp_o_file = temp_file(".tsv");
	file_put_contents($temp_o_file,$output_lines);
	$output = temp_file(".tsv");
	
	//sort file
	exec2("( head -n1 {$temp_o_file}; tail -n+2 {$temp_o_file} | sort -k1.4,1 -k2,2 -n )   > {$output} ");
	return $output;
}

//Checks if column exists in GSvar file and removes it
function remove_col($col_name,&$matrix)
{
	$col_indices = array();
	foreach($matrix->getHeaders() as $index => $name)
	{
		if($col_name == $name)
		{
			$col_indices[] = $index;
		}
	}
	
	if(empty($col_indices)) return;
	
	for($i=count($col_indices)-1;$i>=0;--$i)
	{
		$matrix->removeCol($col_indices[$i]);
	}
}

//Returns whether a certain gene is an oncogene or tsg gene according NCG6.0
function ncg_gene_statements($gene)
{
	$result = array("is_oncogene" => "na", "is_tsg" => "na");
	
	$ncg_file = get_path("data_folder") . "/dbs/NCG6.0/NCG6.0_oncogene.tsv";
	
	if(!file_exists($ncg_file))
	{
		trigger_error("Could not find file with NCG6.0 data",E_USER_WARNING);
		return $result;
	}
	
	$handle = fopen($ncg_file,"r");
	while(!feof($handle))
	{
		$line = fgets($handle);
		if(empty($line)) continue;
		if(starts_with($line,"entrez")) continue;
		
		list($entrez,$ncg_gene,$cgc,$vogelstein,$is_oncogene,$is_tsg) = explode("\t",trim($line));
		
		if(($is_oncogene) == "")
		{
			echo "empty: $entrez $is_oncogene \n";
			die();
		}
		
		if(trim($gene) == trim($ncg_gene))
		{
			$result["is_oncogene"] = $is_oncogene;
			
			$result["is_tsg"] = $is_tsg;
			break;
		}
	}
	
	fclose($handle);
	
	return $result;
}


/********
 * MAIN *
 ********/

if(isset($cgi_snv_in))
{
	$gsvar_input = Matrix::fromTSV($gsvar_in);
	/******************************
	 * REMOVE OLD CGI ANNOTATIONS *
	 ******************************/ 
	//old annotations created by former script
	$old_annotations = ["CGI_drug_assoc","CGI_evid_level","CGI_transcript","CGI_id","CGI_driver_statement","CGI_gene_role","CGI_gene","CGI_consequence"];
	foreach($old_annotations as $annotation_name) remove_col($annotation_name,$gsvar_input);

	/***************************
	 * REMOVE OLD CGI COMMENTS *
	 ***************************/
	function remove_comment($old_comment,&$matrix)
	{
		foreach($matrix->getComments() as $comment)
		{
			if(strpos($comment,$old_comment) !== false) $matrix->removeComment($comment);
		}
	}
	$old_comments = ["CGI_ICD10_CODE","CGI_ICD10_TEXT","CGI_ICD10_DIAGNOSES","GENES_FOR_REIMBURSEMENT","CGI_CANCER_TYPE","CGI_id","CGI_driver_statement","CGI_gene_role","CGI_transcript","CGI_gene","CGI_consequence"];
	foreach($old_comments as $old_comment) remove_comment($old_comment,$gsvar_input);

	//indices of input gsvar position columns
	$i_gsvar_snvs_chr = $gsvar_input->getColumnIndex("chr");
	$i_gsvar_snvs_start = $gsvar_input->getColumnIndex("start");
	$i_gsvar_snvs_genes = $gsvar_input->getColumnIndex("gene");

	//VCF file with CGI mutation analysis
	$cgi_snvs = file(cgi_variant_tsv_to_sorted_vcf($cgi_snv_in),FILE_IGNORE_NEW_LINES);

	//arrays for CGI data
	$col_cgi_id = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_driver_statement = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_gene_role = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_transcript = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_gene = array_fill(0,$gsvar_input->rows(),"");
	$col_cgi_consequence = array_fill(0,$gsvar_input->rows(),"");

	//annotate GSvar file with CGI data, use positions of SNV for annotation
	for($i=0;$i<$gsvar_input->rows();$i++)
	{
		$gsvar_chr = chr_trim($gsvar_input->get($i,$i_gsvar_snvs_chr));
		$gsvar_pos = $gsvar_input->get($i,$i_gsvar_snvs_start);
		
		foreach($cgi_snvs as $cgi_line)
		{
			
			if(starts_with($cgi_line,"#")) continue;
			list($cgi_chr,$cgi_pos,$cgi_ref,$cgi_alt,$cgi_statements) = explode("\t",$cgi_line);
			$cgi_statements = explode(",",$cgi_statements);
			$cgi_chr = chr_trim($cgi_chr);
			
			if((is_numeric($cgi_chr) && is_numeric($gsvar_chr)) && ($cgi_chr > $gsvar_chr)) break;
			if($cgi_chr != $gsvar_chr) continue;
			if(strlen($cgi_ref) == 1 && $gsvar_pos != $cgi_pos) continue;
			if(strlen($cgi_ref) > 1 && $gsvar_pos != $cgi_pos + 1) continue;
			if(strlen($cgi_ref) == 1) //SNPs
			{
				if($gsvar_chr == $cgi_chr and $gsvar_pos == $cgi_pos)
				{
					$col_cgi_id[$i] = $cgi_statements[0];
					$col_cgi_driver_statement[$i] = $cgi_statements[1];
					$col_cgi_gene_role[$i] = $cgi_statements[2];
					$col_cgi_transcript[$i] = $cgi_statements[3];
					$col_cgi_gene[$i] = $cgi_statements[4];
					$col_cgi_consequence[$i] = $cgi_statements[5];
					break;
				}
			}
			elseif(strlen($cgi_ref) > 1) //indels
			{
				//VCF Ref/Alt GXXXXX/G
				//GSVar Ref/Alt XXXXX/-
				$cgi_start_indel = $cgi_pos+1;
				$cgi_ref_indel = substr($cgi_ref,1);
				if($gsvar_chr == $cgi_chr and $gsvar_pos == $cgi_start_indel)
				{
					$col_cgi_id[$i] = $cgi_statements[0];
					$col_cgi_driver_statement[$i] = $cgi_statements[1];
					$col_cgi_gene_role[$i] = $cgi_statements[2];
					$col_cgi_transcript[$i] = $cgi_statements[3];
					$col_cgi_gene[$i] = $cgi_statements[4];
					$col_cgi_consequence[$i] = $cgi_statements[5];
					break;
				}
			}
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

	$gsvar_input->toTSV($out);
}

if($include_ncg)
{
	$gsvar_input = Matrix::fromTSV($gsvar_in);
	

	//Remove potential old NCG annotation
	remove_col("ncg_oncogene",$gsvar_input);
	remove_col("ncg_tsg",$gsvar_input);
	$gsvar_input->toTSV($out);
	
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
	
	$gsvar_input->toTSV($out);
}

?>
