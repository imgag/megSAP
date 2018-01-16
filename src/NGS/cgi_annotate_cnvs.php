<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("cgi_annotate_cnvs", "Annotates CGI CNV data to cnv.tsv file.");
$parser->addInfile("cnv_in", "Input .gsvar-file with SNV data.", false);
$parser->addInfile("cnv_in_cgi", "Input CGI data with SNV annotations",false);
$parser->addOutfile("out", "Output file name", false);
extract($parser->parse($argv));
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

$cnv_region_copy_numbers = $cnv_input->getCol($cnv_input->getColumnIndex("region_copy_numbers"));

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
	$cgi_gene_role = $cgi_gene_roles[$i];
	$cgi_alteration_type = $cgi_alteration_types[$i];
	
	for($j=0;$j<count($cnv_genes);$j++)
	{
		$genes = explode(',',$cnv_genes[$j]);
		
		//calc Copy number alteration type (amp or del) (!per region!) in input cnv file		
		$copy_numbers_region = explode(',',$cnv_region_copy_numbers[$j]);
		$median_copy_number = median($copy_numbers_region);
		$cnv_alteration_type;
		if($median_copy_number>2.)
		{
			$cnv_alteration_type = "AMP";
		} else {
			$cnv_alteration_type = "DEL";
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

$output = $cnv_input;

//remove old CGI data

if($output->getColumnIndex("CGI_drug_assoc",false,false) !== false)
{
	$output->removeCol($output->getColumnIndex("CGI_drug_assoc"));
}
if($output->getColumnIndex("CGI_evid_level",false,false) !== false)
{
	$output->removeCol($output->getColumnIndex("CGI_evid_level"));
}
if($output->getColumnIndex("CGI_genes",false,false) !== false)
{
	$output->removeCol($output->getColumnIndex("CGI_genes"));
}
if($output->getColumnIndex("CGI_driver_statement",false,false) !== false)
{
	$output->removeCol($output->getColumnIndex("CGI_driver_statement"));
}
if($output->getColumnIndex("CGI_gene_role",false,false) !== false)
{
	$output->removeCol($output->getColumnIndex("CGI_gene_role"));
}

//add CGI data
$cancertype = $cgi_input->get(0,$cgi_input->getColumnIndex("cancer"));
$output->addCol($new_genes,"CGI_genes","Genes which were included in CancerGenomeInterpreter.org analysis.");
$output->addCol($new_driver_statement,"CGI_driver_statement","Driver statement for cancer type $cancertype according CancerGenomeInterpreter.org");
$output->addCol($new_gene_role,"CGI_gene_role","Gene role for cancer type $cancertype according CancerGenomeInterpreter.org");
$output->toTSV($out);

?>
