<?php
/** 
	@page db_delete_sample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_delete_sample", "Deletes processed samples data from file system and NGSD.");
$parser->addStringArray("samples", "Processed sample name list.", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db_conn = DB::getInstance($db);
$data_folder = get_path("data_folder");

//remove 
foreach($samples as $ps_name)
{
	$ps_id = get_processed_sample_id($db_conn, $ps_name, false);		
	if ($ps_id==-1)
	{
		print "{$ps_name}: not found in NGSD!\n";
	}
	else
	{
		print "{$ps_name}: deleting...\n";

		
		//file system: remove sample folder
		$ps_info = get_processed_sample_info($db_conn, $ps_name);
		$ps_folder = $ps_info['ps_folder'];
		print "    folder $ps_folder\n";
		exec2("rm -rf {$ps_folder}");
		
		//file system: remove coverage data
		$cov_file = "{$data_folder}/coverage/".$ps_info['sys_name_short']."/{$ps_name}.cov";
		print "    coverage data $cov_file\n";
		exec2("rm -rf {$cov_file}");
		
		print "    data from NGSD\n";
		
		//NGSD: remove FK tables
		$tables = array("detected_variant"=>"processed_sample_id", "detected_somatic_variant"=>"processed_sample_id_tumor", "detected_somatic_variant"=>"processed_sample_id_normal", "processed_sample_qc"=>"processed_sample_id", "diag_status"=>"processed_sample_id", "kasp_status"=>"processed_sample_id");
		foreach($tables as $table => $col)
		{
			$db_conn->executeStmt("DELETE FROM {$table} WHERE {$col}='$ps_id'");
		}
		
		//NGSD: remove analysis jobs
		$job_ids = $db_conn->getValues("SELECT j.id FROM analysis_job j, analysis_job_sample js WHERE js.analysis_job_id=j.id AND js.processed_sample_id='$ps_id'");
		foreach($job_ids as $job_id)
		{
			$db_conn->executeStmt("DELETE FROM analysis_job_sample WHERE analysis_job_id='$job_id'");
			$db_conn->executeStmt("DELETE FROM analysis_job_history WHERE analysis_job_id='$job_id'");
			$db_conn->executeStmt("DELETE FROM analysis_job WHERE id='$job_id'");
		}
		
		//NGSD: remove processed sample
		$db_conn->executeStmt("DELETE FROM processed_sample WHERE id='$ps_id'");
	}	
}

?>