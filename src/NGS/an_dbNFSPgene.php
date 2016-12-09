<?php 
/** 
	@page an_dbNFSPgene	
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_dbNFSPgene", "\$Rev: 745 $", "Annotate additional columns from dbNFSPgene.");
$parser->addInfile("in",  "Input file in tsv format.", false);
$parser->addOutfile("out", "Output file in tsv format.", false);
//optional
$parser->addInt("idx_coding", "Coding column (-1 = auto select - requires column 'coding_and_splicing').", true, -1);
extract($parser->parse($argv));

//read dbNFSPgene and generate map of refseq-IDs
$db = get_path("data_folder")."/dbs/dbNSFP/dbNSFP2.9_gene";
$refseq_map = array();
$dbNFSPgene = Matrix::fromTSV($db);
$dbNFSPgene->setHeaders($dbNFSPgene->getRow(0));
$idx_k = $dbNFSPgene->getColumnIndex("Pathway(KEGG)_full");
$idx_f = $dbNFSPgene->getColumnIndex("Function_description");
$idx_b = $dbNFSPgene->getColumnIndex("Pathway(BioCarta)_full");
for($i=0;$i<$dbNFSPgene->rows();++$i)
{
	$row = $dbNFSPgene->getRow($i);
	
	foreach(explode(",", $row[9]) as $refseq_id)
	{
		if($refseq_id == ".")	continue;
//		if(isset($refseq_map[$refseq_id]))	trigger_error("Found Refseq-ID '$refseq_id' multiple times.",E_USER_NOTICE);
		if(!isset($refseq_map[$refseq_id]))	$refseq_map[$refseq_id] = array();
		$refseq_map[$refseq_id][] = $i;
	}
}

//write info fields; get refseq IDs from vcf-file and add dbNFSPgene information
$variant_file = Matrix::fromTSV($in);
$new_cols = array("KEGG"=>array(),"BioCarta"=>array(),"function"=>array());

if($idx_coding == -1)
{
	$idx_coding = $variant_file->getColumnIndex("coding_and_splicing");
	if($idx_coding == -1)	trigger_error("Could not identify column 'coding_and_splicing'.",E_USER_ERROR);
}


//remove columns from previous annotation
$col = $variant_file->getColumnIndex("Pathway_KEGG_full",false,false);
if($col!==FALSE)	$variant_file->removeCol($col);
$col = $variant_file->getColumnIndex("Function_description",false,false);
if($col!==FALSE)	$variant_file->removeCol($col);
$col = $variant_file->getColumnIndex("Pathway_BioCarta_full",false,false);
if($col!==FALSE)	$variant_file->removeCol($col);

$genes_not_found = array();
for($i=0;$i<$variant_file->rows();++$i)
{
	$row = $variant_file->getRow($i);
	$coding = explode(",",$row[$idx_coding]);
	$kegg = array();
	$biocarta = array();
	$func = array();
	foreach($coding as $c)
	{
		$k = "";
		$f = "";
		$b = "";
		if(!empty($c))
		{
			$parts = explode(":",$c);
			$refseq_id = $parts[1];
			if(strpos($refseq_id, ".")!==FALSE)	$refseq_id = substr($refseq_id,0,strpos($refseq_id,'.'));
			if(empty($refseq_id))	continue;
			if(isset($refseq_map[$refseq_id]))
			{
				foreach($refseq_map[$refseq_id] as $r);
				{
					$k .= $dbNFSPgene->get($r,$idx_k);
					$f .= $dbNFSPgene->get($r,$idx_f);
					$b .= $dbNFSPgene->get($r,$idx_b);
				}
			}
			else
			{
				$genes_not_found[] = $refseq_id;
			}
		}
		$kegg[] = $k;
		$biocarta[] = $b;	
		$func[] = $f;	
	}
	$new_cols["KEGG"][] = trim(implode(",",array_values(array_unique($kegg))),",.");	//17
	$new_cols["function"][] = trim(implode(",",array_values(array_unique($func))),",.");	//18
	$new_cols["BioCarta"][] = trim(implode(",",array_values(array_unique($biocarta))),",.");	//14
}
$genes_not_found = array_unique($genes_not_found);
if(!empty($genes_not_found))	trigger_error("Could not find the following genes: ".implode(", ", $genes_not_found).".",E_USER_WARNING);

//annotate at end if interpro column is not available
$idx_interpro = $variant_file->getColumnIndex("interpro",false,false);
if($idx_interpro===FALSE)	$idx_interpro = $variant_file->cols();
$variant_file->insertCol($idx_interpro+1, $new_cols["KEGG"], "Pathway_KEGG_full");
$variant_file->insertCol($idx_interpro+2, $new_cols["function"], "Function_description");
$variant_file->insertCol($idx_interpro+3, $new_cols["BioCarta"], "Pathway_BioCarta_full");
$variant_file->toTSV($out);

?>