<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("cgi_snvs_to_gsvar", "Annotates CGI data to GSvar file.");
$parser->addInfile("gsvar_in", "Input .gsvar-file with SNV data.", false);
$parser->addInfile("cgi_snv_in", "Input CGI data with SNV annotations",false);
$parser->addOutfile("out", "Output file name", false);
extract($parser->parse($argv));
$gsvar_input = Matrix::fromTSV($gsvar_in);


//remove old CGI annotations if annotated 

if($gsvar_input->getColumnIndex("CGI_drug_assoc",false,false) !== false)
{
	$gsvar_input->removeCol($gsvar_input->getColumnIndex("CGI_drug_assoc"));
	$gsvar_input->removeCol($gsvar_input->getColumnIndex("CGI_evid_level"));
	$gsvar_input->removeCol($gsvar_input->getColumnIndex("CGI_transcript"));
}
//remove old columns created by this script
if($gsvar_input->getColumnIndex("CGI_driver_statement",false,false) !== false)
{
	$gsvar_input->removeCol($gsvar_input->getColumnIndex("CGI_driver_statement",false,false));
}
if($gsvar_input->getColumnIndex("CGI_gene_role",false,false) !== false)
{
	$gsvar_input->removeCol($gsvar_input->getColumnIndex("CGI_gene_role",false,false));
}
if($gsvar_input->getColumnIndex("CGI_transcript",false,false) !== false)
{
	$gsvar_input->removeCol($gsvar_input->getColumnIndex("CGI_transcript",false,false));
}
if($gsvar_input->getColumnIndex("CGI_gene",false,false) !== false)
{
	$gsvar_input->removeCol($gsvar_input->getColumnIndex("CGI_gene",false,false));
}
if($gsvar_input->getColumnIndex("CGI_consequence",false,false) !== false)
{
	$gsvar_input->removeCol($gsvar_input->getColumnIndex("CGI_consequence",false,false));
}

//indices of gsvar position columns
$i_gsvar_snvs_chr = $gsvar_input->getColumnIndex("chr");
$i_gsvar_snvs_start = $gsvar_input->getColumnIndex("start");
$i_gsvar_snvs_genes = $gsvar_input->getColumnIndex("gene");


//tsv with cgi results
$cgi_snvs = Matrix::fromTSV($cgi_snv_in);


//column indices of input snv file
$i_cgi_snvs_pos = $cgi_snvs->getColumnIndex("input");
$i_cgi_snvs_driver = $cgi_snvs->getColumnIndex("driver");
$i_cgi_snvs_driver_statement = $cgi_snvs->getColumnIndex("driver_statement");
$i_cgi_snvs_gene_role = $cgi_snvs->GetColumnIndex("gene_role");
$i_cgi_snvs_genes = $cgi_snvs->GetColumnIndex("gene");
$i_cgi_snvs_transcript = $cgi_snvs->GetColumnIndex("transcript");
$i_cgi_snvs_consequence = $cgi_snvs->GetColumnIndex("consequence");

//CGI columns for positions
$cgi_snvs_start = array();
$cgi_snvs_chr = array();
//reference allele in vcf file, important to annotate indels
$cgi_snvs_ref_allele = array();
//cgi gene names
$cgi_snvs_gene_names = array();
$cgi_snvs_transcript = array();
$cgi_snvs_consequence = array();

//fill CGI columns for important results
for($i=0;$i<$cgi_snvs->rows();$i++)
{
	$coordinates = explode('|',$cgi_snvs->getCol($i_cgi_snvs_pos)[$i]);
	$cgi_snvs_chr[] = "chr".$coordinates[0];
	$cgi_snvs_start[] = $coordinates[1];
	$cgi_snvs_ref_allele[] = $coordinates[2];
}
$cgi_snvs_gene_names = $cgi_snvs->getCol($i_cgi_snvs_genes);
//update gene names in CGI file
$cgi_snvs_gene_names = approve_gene_names($cgi_snvs_gene_names);

$cgi_snvs_transcript = $cgi_snvs->getCol($i_cgi_snvs_transcript);
$cgi_snvs_consequence = $cgi_snvs->getCol($i_cgi_snvs_consequence);

//new columns in GSVar file
$gsvar_snvs_new_driver = array();
$gsvar_snvs_new_driver_statement = array();
$gsvar_snvs_gene_role = array();
$gsvar_snvs_transcript = array();
$gsvar_snvs_new_cgi_gene = array();
$gsvar_snvs_cgi_consequence = array();

for($i=0;$i<$gsvar_input->rows();$i++)
{
	$gsvar_snvs_new_driver[] = "";
	$gsvar_snvs_new_driver_statement[] = "";
	$gsvar_snvs_gene_role[] = "";
	$gsvar_snvs_transcript[] = "";
	$gsvar_snvs_new_cgi_gene[] = "";
	$gsvar_snvs_cgi_consequence[] = "";
}



//annotate GSvar file with CGI data, use positions of SNV for annotation
for($i=0;$i<$gsvar_input->rows();$i++)
{
	$gsvar_chr = $gsvar_input->get($i,$i_gsvar_snvs_chr);
	$gsvar_start = $gsvar_input->get($i,$i_gsvar_snvs_start);
	$gsvar_gene = $gsvar_input->get($i,$i_gsvar_snvs_genes);
	
	for($j=0;$j<$cgi_snvs->rows();$j++)
	{
		$cgi_chr = $cgi_snvs_chr[$j];
		$cgi_start = $cgi_snvs_start[$j];
		$cgi_ref = $cgi_snvs_ref_allele[$j];
		
		
		$cgi_gene = $cgi_snvs_gene_names[$j];
		$cgi_transcript = $cgi_snvs_transcript[$j];
		
		$cgi_snvs_driver = $cgi_snvs->getCol($i_cgi_snvs_driver)[$j];
		$cgi_snvs_driver_statement = $cgi_snvs->getCol($i_cgi_snvs_driver_statement)[$j];
		$cgi_snvs_gene_role = $cgi_snvs->getCol($i_cgi_snvs_gene_role)[$j];
		$cgi_consequence = $cgi_snvs_consequence[$j];
		
		//Distinguish between SNVs and indels
		if(strlen($cgi_ref) == 1)
		{
			if($gsvar_chr == $cgi_chr and $gsvar_start == $cgi_start)
			{
				//check whether gene names from CGI and in GSVAR match
				if(strpos($gsvar_gene,$cgi_gene) === false)
				{
					trigger_error("gene in CGI does not match GSVar gene: $cgi_gene != $gsvar_gene",E_USER_WARNING);
				}
				
				$gsvar_snvs_new_driver[$i] = $cgi_snvs_driver;
				$gsvar_snvs_new_driver_statement[$i] = $cgi_snvs_driver_statement;
				$gsvar_snvs_gene_role[$i] = $cgi_snvs_gene_role;
				$gsvar_snvs_transcript[$i] = $cgi_transcript;
				$gsvar_snvs_new_cgi_gene[$i]= $cgi_gene;
				$gsvar_snvs_cgi_consequence[$i] = $cgi_consequence;
				break;
			}
		}
		elseif(strlen($cgi_ref)> 1) //indels, format in vcf is different than in GSVar
		{
			//VCF Ref/Alt GXXXXX/G
			//GSVar Ref/Alt XXXXX/-
			$cgi_start_indel = $cgi_start+1;
			$cgi_ref_indel = substr($cgi_ref,1);
			if($gsvar_chr == $cgi_chr and $gsvar_start == $cgi_start_indel)
			{
				$gsvar_snvs_new_driver[$i] = $cgi_snvs_driver;
				$gsvar_snvs_new_driver_statement[$i] = $cgi_snvs_driver_statement;
				$gsvar_snvs_gene_role[$i] = $cgi_snvs_gene_role;
				$gsvar_snvs_transcript[$i] = $cgi_transcript;
				$gsvar_snvs_new_cgi_gene[$i] = $cgi_gene;
				$gsvar_snvs_cgi_consequence[$i] = $cgi_consequence;
				break;
			}
		}
	}
}

//get cancer type from CGI input file
$cancer_type_cgi = $cgi_snvs->get(0,$cgi_snvs->getColumnIndex("cancer"));

//insert header line which describes genes neccessary for reimbursement with health insurance
$dictionary = Matrix::fromTSV(repository_basedir()."/data/dbs/Ontologies/icd10_cgi_dictionary.tsv");
$row_cancer_type = -1;
$i_cancer_acronym = $dictionary->getColumnIndex("cgi_acronym");
for($i=0;$i<$dictionary->rows();$i++)
{
	if($dictionary->get($i,$i_cancer_acronym) == $cancer_type_cgi)
	{
		$row_cancer_type = $i;
	}
}

$icd10_diagnosis_code = "";
$icd10_diagnosis_text = "";
if($row_cancer_type == -1)
{
	trigger_error("Could not determine ICD-10 diagnosis.",E_USER_WARNING);
}
else
{
	$icd10_diagnosis_code = $dictionary->get($row_cancer_type,$dictionary->getColumnIndex("icd10_code"));
	$icd10_diagnosis_text = ($dictionary->get($row_cancer_type,$dictionary->getColumnIndex("desc_de")));
}

$output_comments = $gsvar_input->getComments();
//add information about Diagnoses to GSVar file

//remove old comments
for($i=0;$i<count($output_comments);$i++)
{
	if(strpos($output_comments[$i],"CGI_ICD10_CODE") !== false)
	{
		$gsvar_input->removeComment($output_comments[$i]);
	}
}
for($i=0;$i<count($output_comments);$i++)
{
	if(strpos($output_comments[$i],"CGI_ICD10_TEXT") !== false)
	{
		$gsvar_input->removeComment($output_comments[$i]);
	}
}
for($i=0;$i<count($output_comments);$i++)
{
	if(strpos($output_comments[$i],"CGI_ICD10_DIAGNOSES") !== false)
	{
		$gsvar_input->removeComment($output_comments[$i]);
	}
}
for($i=0;$i<count($output_comments);$i++)
{
	if(strpos($output_comments[$i],"GENES_FOR_REIMBURSEMENT") !== false)
	{
		$gsvar_input->removeComment($output_comments[$i]);
	}
}
for($i=0;$i<count($output_comments);$i++)
{
	if(strpos($output_comments[$i],"CGI_CANCER_TYPE") !== false)
	{
		$gsvar_input->removeComment($output_comments[$i]);
	}
}
$gsvar_input->addComment("#CGI_ICD10_TEXT=$icd10_diagnosis_text");
$gsvar_input->addComment("#CGI_ICD10_CODE=$icd10_diagnosis_code");
$gsvar_input->addComment("#CGI_CANCER_TYPE=$cancer_type_cgi");

$gsvar_input->addCol($gsvar_snvs_new_driver_statement,"CGI_driver_statement","Oncogenic Classification according CGI for tumor type $cancer_type_cgi");
$gsvar_input->addCol($gsvar_snvs_gene_role,"CGI_gene_role","CGI gene role. LoF: Loss of Function, Act: Activating");
$gsvar_input->addCol($gsvar_snvs_transcript,"CGI_transcript","CGI Ensembl transcript ID");
$gsvar_input->addCol($gsvar_snvs_new_cgi_gene,"CGI_gene","Gene returned by CGI");
$gsvar_input->addCol($gsvar_snvs_cgi_consequence,"CGI_consequence","Consequence of the mutation assessed by CGI");
$gsvar_input->toTSV($out);
?>
