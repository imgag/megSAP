<?php
/** 
	@page export_ngsd 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function export_table(&$out_handle, $table, $conditions = "1")
{
	//init
	global $db;
	
	print "  Exporting table {$table}";
	
	//extract field names
	$fields = array();
	$res = $db->executeQuery("DESCRIBE {$table}");
	foreach($res as $row)
	{
		$field = $row['Field'];
		$fields[$field] = array(
							'null'=> ($row['Null']=="YES")
							);
	}
	
	//write output
	$max_elements = 1000;
	$res = $db->executeQuery("SELECT * FROM {$table} WHERE {$conditions}");
	$row_count = count($res);
	if ($row_count>0)
	{
		for ($i=0; $i<$row_count; ++$i)
		{
			if ($i%$max_elements==0)
			{
				fputs($out_handle, "\n");
				fputs($out_handle, "INSERT INTO `{$table}`(`".implode("`, `", array_keys($fields))."`) VALUES\n");
			}
			$row = $res[$i];
			$line = "  (";
			
			$items = array();
			foreach($fields as $field => $field_meta)
			{
				$item = $row[$field];
				if ($item=="" && $field_meta['null'])
				{
					$items[] = "NULL";
				}
				else
				{
					$items[] = $db->quote($item);
				}
			}
			$line .= implode(",", $items);
			$line .= ")";
			$line .= ($i+1==$row_count || ($i+1)%$max_elements==0) ? ";" : ",";
			fputs($out_handle, $line."\n");
		}
	}
	
	print " ({$row_count} rows)\n";
}

//parse command line arguments
$parser = new ToolBase("export_ngsd", "Exports information from NGSD.");
$parser->addOutfile("out", "Output file that SQL queries are written to.", false);
$parser->addStringArray("samples", "Sample names to include in the export.", true);
$parser->addEnum("db",  "Database to export from to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db = DB::getInstance($db);
if (file_exists($out)) unlink($out);
$out_handle = fopen($out, "a");

//remove all initial data from NGSD (not by TUNCATE in case someone applies the SQL import to a production database)
fputs($out_handle, "DELETE FROM user WHERE user_id='admin' OR user_id='genlab_import';\n");
fputs($out_handle, "DELETE FROM species WHERE name='human';\n");
fputs($out_handle, "DELETE FROM genome WHERE build='hg19' OR build='hg38';\n");

//export base tables we always need
print "Exporting complete base tables...\n";
export_table($out_handle, "user");
export_table($out_handle, "device");
export_table($out_handle, "gene");
export_table($out_handle, "gene_alias");
export_table($out_handle, "gene_transcript");
export_table($out_handle, "gene_exon");
export_table($out_handle, "geneinfo_germline");
export_table($out_handle, "genome");
export_table($out_handle, "hpo_term");
export_table($out_handle, "hpo_parent");
export_table($out_handle, "hpo_genes");
export_table($out_handle, "mid");
export_table($out_handle, "omim_gene");
export_table($out_handle, "omim_phenotype");
export_table($out_handle, "processing_system");
export_table($out_handle, "project");
export_table($out_handle, "qc_terms");
export_table($out_handle, "sender");
export_table($out_handle, "sequencing_run");
export_table($out_handle, "runqc_read");
export_table($out_handle, "runqc_lane");
export_table($out_handle, "species");

//export sample data
foreach($samples as $sample)
{
	print "\n";
	print "Exporting sample {$sample}...\n";
	
	//sample
	$sample_id = $db->getValue("SELECT id FROM sample WHERE name='{$sample}'");
	export_table($out_handle, "sample", "id={$sample_id}");
	
	//sample disease information
	export_table($out_handle, "sample_disease_info", "sample_id={$sample_id}");
	
	//processed samples
	$processed_sample_ids = $db->getValues("SELECT id FROM processed_sample WHERE sample_id={$sample_id}");
	export_table($out_handle, "processed_sample", "id IN (".implode(",", $processed_sample_ids).")");
	
	//processed sample QC
	export_table($out_handle, "processed_sample_qc", "processed_sample_id IN (".implode(",", $processed_sample_ids).")");
	
	//diagnostic status
	export_table($out_handle, "diag_status", "processed_sample_id IN (".implode(",", $processed_sample_ids).")");
	
	//cnv callset
	export_table($out_handle, "cnv_callset", "processed_sample_id IN (".implode(",", $processed_sample_ids).")");
	$cnv_callset_ids = $db->getValues("SELECT id FROM cnv_callset WHERE processed_sample_id IN (".implode(",", $processed_sample_ids).")");
	export_table($out_handle, "cnv", "cnv_callset_id IN (".implode(",", $cnv_callset_ids).")");
}

fclose($out_handle);

?>