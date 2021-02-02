<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("2019_09_diag_status_update", "NGSD update for diag_status: Removes columns genes_causal, inheritance_mode and genes_incidental; stores the date in comments.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$marker = "Diagnostik status infos (old):";
$db = DB::getInstance($db);
$user_id = $db->getValue("SELECT id FROM user WHERE user_id='unknown'");

//process diag status
$result = $db->executeQuery("SELECT s.name, ps.process_id, ds.* FROM diag_status ds, sample s, processed_sample ps WHERE ps.sample_id=s.id AND ps.id=ds.processed_sample_id");
foreach($result as $row)
{
	$ps = $row['name']."_".$row['process_id'];
	$ps_id = $row['processed_sample_id'];
	
	$genes_causal = trim($row['genes_causal']);
	if ($genes_causal=="n/a") $genes_causal = "";
	$genes_incidental = trim($row['genes_incidental']);
	if ($genes_incidental=="n/a") $genes_incidental = "";
	$inheritance_mode = trim($row['inheritance_mode']);
	if ($inheritance_mode=="n/a") $inheritance_mode = "";
			

	if ($genes_causal.$genes_incidental.$inheritance_mode!="")
	{
		//update NGSD entry
		$comment = trim($row['comment']);
		if (contains($comment, $marker))
		{
			print "$ps ($ps_id): skipping - already updated\n";
		}
		else
		{
			print "$ps ($ps_id): updating...\n";
			
			if ($comment!="") $comment .= "\n";
			$comment .= "{$marker}\n";
			if ($genes_causal!="") $comment .= "genes_causal: {$genes_causal}\n";
			if ($genes_incidental!="") $comment .= "genes_incidental: {$genes_incidental}\n";
			if ($inheritance_mode!="") $comment .= "inheritance_mode: {$inheritance_mode}\n";
			
			$db->executeStmt("UPDATE diag_status SET comment=:0 WHERE processed_sample_id=:1", array("0"=>$comment, "1"=>$ps_id));
		}
		
		//create report configuration (if possible)
		if ($genes_causal!="")
		{
			$report_id = $db->getValue("SELECT id FROM report_configuration WHERE processed_sample_id='{$ps_id}'", -1);
			if ($report_id>=0)
			{
				print "  report config: skipping - already present\n";
			}
			else
			{
				//get approved gene names
				$tmp = strtr($genes_causal, ";,=/:()-", "        ");
				$tmp = explode(" ", $tmp);
				$genes = array();
				foreach($tmp as $gene)
				{
					$gene = trim($gene);
					if ($gene=="") continue;
					
					
					list($tmp2) = exec2("echo '$gene' | ".get_path("ngs-bits")."GenesToApproved");
					if (contains($tmp2[0], "ERROR:")) continue;
					list($gene_approved) = explode("\t", $tmp2[0]);
					
					$genes[] = $gene_approved;
				}
				
				//get class 4/5 variants in genes
				$variant_ids = array();
				$result2 = $db->executeQuery("SELECT v.id, v.coding FROM variant v, detected_variant dv, variant_classification vc WHERE dv.variant_id=v.id AND dv.processed_sample_id={$ps_id} AND vc.variant_id=v.id AND (vc.class='4' OR vc.class='5')");
				foreach($result2 as $row2)
				{
					foreach($genes as $gene)
					{
						if (contains($row2['coding'], "{$gene}:")) $variant_ids[] = $row2['id'];
					}
				}
				$variant_ids = array_unique($variant_ids);
				
				if (count($variant_ids)>0)
				{
					print "  report config: updating...\n";
					$db->executeStmt("INSERT INTO `report_configuration`(`processed_sample_id`, `created_by`) VALUES ({$ps_id},{$user_id})");
					$report_id = $db->lastInsertId();
					foreach($variant_ids as $v_id)
					{
						$db->executeStmt("INSERT INTO `report_configuration_variant`(`report_configuration_id`, `variant_id`, `type`, `causal`, `inheritance`, `de_novo`, `mosaic`, `compound_heterozygous`, `exclude_artefact`, `exclude_frequency`, `exclude_phenotype`, `exclude_mechanism`, `exclude_other`, `comments`, `comments2`) VALUES ('{$report_id}', '{$v_id}', 'diagnostic variant', TRUE, 'n/a', FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, '', '')");	
					}
				}
			}
		}
	}
	else
	{
		print "$ps ($ps_id): skipped - no info to copy\n";
	}
}

//delete columns that are not used anymore
$db->executeStmt("ALTER TABLE `diag_status` DROP `genes_causal`, DROP `inheritance_mode`, DROP `genes_incidental`;");


?>
