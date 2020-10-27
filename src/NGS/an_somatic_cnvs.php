<?php
/** 
	@page an_somatic_cnvs
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_somatic_cnvs", "Annotates additional somatic data to ClinCNV file (CGI / NCG6.0 / RNA counts).");
$parser->addInfile("cnv_in", "Input CNV file.", false);
$parser->addInfile("cnv_in_cgi", "Input CNV CGI data",true);
$parser->addFlag("include_ncg", "Annotate column with info from NCG6.0 whether a gene is TSG or oncogene");
$parser->addFlag("include_cytoband", "Annotate column with cytoband names.");
$parser->addInfile("rna_counts", "Annotate transcript counts from RNA counts file", true);
$parser->addString("rna_id", "Processed sample ID of RNA to be annotated.", true);
$parser->addString("rna_ref_tissue", "RNA tissue type for annotation of RNA reference transciptome data.", true);
$parser->addOutfile("out", "Output file name", false);
extract($parser->parse($argv));

if(!isset($cnv_in_cgi) && !$include_ncg && !isset($rna_counts) && !$include_cytoband)
{
	trigger_error("At least one input file with annotation data must be given.", E_USER_ERROR);
}

$cnv_input = Matrix::fromTSV($cnv_in);


if(isset($cnv_in_cgi))
{
	$cnv_input = Matrix::fromTSV($cnv_in);
	$cgi_input = Matrix::fromTSV($cnv_in_cgi);
	
	//data from CGI cnv file
	$cgi_genes = $cgi_input->getCol($cgi_input->getColumnIndex("gene"));
	$cgi_driver_statements = $cgi_input->getCol($cgi_input->getColumnIndex("driver_statement"));
	$cgi_gene_roles = $cgi_input->getCol($cgi_input->getColumnIndex("gene_role"));
	$cgi_alteration_types = $cgi_input->getCol($cgi_input->getColumnIndex("cna"));
	//data from cnv file
	$cnv_genes = $cnv_input->getCol($cnv_input->getColumnIndex("genes"));
	//approved cnv genes
	$approved_cnv_genes = array();
	foreach($cnv_genes as $gene)
	{
		$temp_genes = approve_gene_names(explode(',',$gene));
		$approved_cnv_genes[] = implode(',',$temp_genes);
	}

	$i_cn_cnvhunter = $cnv_input->getColumnIndex("region_copy_numbers",false,false); //CNVHunter
	$i_cn_clincnv = $cnv_input->getColumnIndex("CN_change",false,false); //ClinCnv

	if(($i_cn_cnvhunter == false && $i_cn_clincnv == false) || ($i_cn_cnvhunter !== false && $i_cn_clincnv !== false))
	{
		trigger_error("Unknown format of CNV file {$cnv_in}. Aborting...",E_USER_ERROR);
	}

	$is_cnvhunter = false;
	if($i_cn_clincnv == false) $is_cnvhunter = true;

	$cn = $cnv_input->getCol($i_cn_clincnv); //column with copy numbers
	if($is_cnvhunter) $cn = $cnv_input->getCol($i_cn_cnvhunter); //CNVhunter file

	//new columns for CNV input
	$new_driver_statement = array();
	$new_gene_role = array();
	$new_genes = array();
	for($i=0;$i<$cnv_input->rows();$i++)
	{
		$new_driver_statement[] = "";
		$new_gene_role[] = "";
		$new_genes[] = "";
	}

	for($i=0;$i<count($cgi_genes);$i++)
	{
		$cgi_gene = $cgi_genes[$i];
		$cgi_driver_statement = $cgi_driver_statements[$i];
		//In rare cases, CGI statement contains ",": must be removed because it is used as separator
		$cgi_driver_statement = str_replace(",",";",$cgi_driver_statement);
		
		$cgi_gene_role = $cgi_gene_roles[$i];
		$cgi_alteration_type = $cgi_alteration_types[$i];
		
		for($j=0;$j<count($cnv_genes);$j++)
		{
			$genes = explode(',',$cnv_genes[$j]);
			
			//calc Copy number alteration type (amp or del) (!per region!) in input cnv file		

			$cnv_alteration_type = "";
			
			if($is_cnvhunter) //CNVHunter
			{
				$copy_numbers_region = explode(',',$cn[$j]);
				$median_copy_number = median($copy_numbers_region);
				if($median_copy_number>2.)
				{
					$cnv_alteration_type = "AMP";
				} else {
					$cnv_alteration_type = "DEL";
				}
			}
			else //ClinCNV
			{
				if($cn[$j] > 2.) $cnv_alteration_type = "AMP";
				elseif($cn[$j] < 2.) $cnv_alteration_type = "DEL";
				else $cnv_alteration_type = "NA";
			}
			
			foreach($genes as $cnv_gene)
			{
				if($cgi_gene == $cnv_gene)
				{	
					//alteration types must match
					if($cgi_alteration_type != $cnv_alteration_type)
					{
						continue;
					}
					
					if($new_genes[$j] == "")
					{
						$new_genes[$j] = $cgi_gene;
						$new_driver_statement[$j] = $cgi_driver_statement;
						$new_gene_role[$j] = $cgi_gene_role;
					} else {
						$new_genes[$j] = $new_genes[$j].','.$cgi_gene;
						$new_driver_statement[$j] = $new_driver_statement[$j] .','. $cgi_driver_statement;
						$new_gene_role[$j] = $new_gene_role[$j] . ',' . $cgi_gene_role;					
					}
				}
			}
			
		}
	}



	//remove old CGI data

	if($cnv_input->getColumnIndex("CGI_drug_assoc",false,false) !== false)
	{
		$cnv_input->removeCol($cnv_input->getColumnIndex("CGI_drug_assoc"));
	}
	if($cnv_input->getColumnIndex("CGI_evid_level",false,false) !== false)
	{
		$cnv_input->removeCol($cnv_input->getColumnIndex("CGI_evid_level"));
	}
	if($cnv_input->getColumnIndex("CGI_genes",false,false) !== false)
	{
		$cnv_input->removeCol($cnv_input->getColumnIndex("CGI_genes"));
	}
	if($cnv_input->getColumnIndex("CGI_driver_statement",false,false) !== false)
	{
		$cnv_input->removeCol($cnv_input->getColumnIndex("CGI_driver_statement"));
	}
	if($cnv_input->getColumnIndex("CGI_gene_role",false,false) !== false)
	{
		$cnv_input->removeCol($cnv_input->getColumnIndex("CGI_gene_role"));
	}

	//add CGI data
	$cancertype = $cgi_input->get(0,$cgi_input->getColumnIndex("cancer"));

	$comments = $cnv_input->getComments();
	for($i=0;$i<count($comments);++$i)
	{
		if(strpos($comments[$i],"#CGI_CANCER_TYPE") !== false)
		{
			$cnv_input->removeComment($comments[$i]);
		}
	}
	$cnv_input->addComment("#CGI_CANCER_TYPE={$cancertype}");

	$cnv_input->addCol($new_genes,"CGI_genes","Genes which were included in CancerGenomeInterpreter.org analysis.");
	$cnv_input->addCol($new_driver_statement,"CGI_driver_statement","Driver statement for cancer type $cancertype according CancerGenomeInterpreter.org");
	$cnv_input->addCol($new_gene_role,"CGI_gene_role","Gene role for cancer type $cancertype according CancerGenomeInterpreter.org");
}

if($include_ncg)
{

	$i_cgi_genes = $cnv_input->getColumnIndex("CGI_genes",false,false);
	
	if($i_cgi_genes === false) 
	{
		trigger_error("Cannot annotate file $cnv_in with TCG6.0 data because there is no column CGI_genes.",E_USER_WARNING);
		exit(1);
	}
	
	$cnv_input->removeColByName("ncg_oncogene");
	$cnv_input->removeColByName("ncg_tsg");
	
	$col_tsg = array();
	$col_oncogene = array();
	
	//Annotate per CGI gene
	for($row=0;$row<$cnv_input->rows();++$row)
	{
		$cgi_genes = explode(",",trim($cnv_input->get($row,$i_cgi_genes)));
		
		$tsg_statements = array();
		$oncogene_statements = array();
		
		foreach($cgi_genes as $cgi_gene)
		{
			$statement = ncg_gene_statements($cgi_gene);
			$tsg_statements[] = $statement["is_tsg"];
			$oncogene_statements[] = $statement["is_oncogene"];
		}
		$col_oncogene[] = implode(",",$oncogene_statements);
		$col_tsg[] = implode(",",$tsg_statements);
	}
	$cnv_input->addCol($col_oncogene,"ncg_oncogene","1:gene is oncogene according NCG6.0, 0:No oncogene according NCG6.0, na: no information available about gene in NCG6.0. Order is the same as in column CGI_genes.");
	$cnv_input->addCol($col_tsg,"ncg_tsg","1:gene is TSG according NCG6.0, 0:No TSG according NCG6.0, na: no information available about gene in NCG6.0. Order is the same as in column CGI_genes.");
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
		$temp = approve_gene_names(explode(",",trim($cnv_input->get($row, $i_genes))));
		foreach($temp as $gene) $genes_of_interest[] = $gene;
	}
	
	//Create result array of genes and tpm that occur in RNA_counts and CNV file
	$i_rna_genes = -1;
	$i_rna_tpm = -1;
	$results = array();
	$handle = fopen($rna_counts, "r");
	while(!feof($handle))
	{
		$line = fgets($handle);
		
		//Determine column indices
		if(starts_with($line,"#gene_id"))
		{
			$parts = explode("\t",$line);
			for($i=0; $i<count($parts); ++$i)
			{
				if($parts[$i] == "gene_name") $i_rna_genes = $i;
				if($parts[$i] == "tpm") $i_rna_tpm = $i;
			}
		}
		
		if(starts_with($line,"#")) continue;
		if(empty($line)) continue;
		
		list(,,,,$tpm, $rna_gene) = explode("\t", $line);
		
		if(in_array($rna_gene, $genes_of_interest))
		{
			$results[$rna_gene] = $tpm;
		}
	}
	fclose($handle);
	
	$ref_file = get_path("data_folder") . "/dbs/gene_expression/rna_tissue_hpa.tsv";
	$ref_results = array();
	$handle = fopen($ref_file,"r");
	
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