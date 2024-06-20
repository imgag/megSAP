<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");
require_once("framework.php");

$name = "merge_samples";

function create_folder($name, $files_to_copy, $dummy_file = "")
{
	$folder_path = data_folder()."/+analysis/merge_samples/".$name."/";
	if (! is_dir($folder_path))
	{
		mkdir($folder_path, 0744, true);
	}
	
	foreach($files_to_copy as $file=>$new_name)
	{
		copy($file, $folder_path.$new_name);
	}
	
	//create empty dummy file
	if ($dummy_file != "")
	{
		file_put_contents($folder_path.$dummy_file, "");
	}
}

function clear_analysis_folder()
{
	if (is_dir(data_folder()."/+analysis/merge_samples/"))
	{
		exec("rm -r ".data_folder()."/+analysis/merge_samples/");
	}
}

function reset_test_db(&$db, $tool_name)
{
	check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$tool_name}.sql");
	
	$base_folder_path = data_folder()."/+analysis/merge_samples/";
	
	for($id=1; $id<17; $id++)
	{
		$sample_name_parts = $db->executeQuery("SELECT s.name, ps.process_id FROM processed_sample as ps LEFT JOIN sample as s ON s.id=ps.sample_id WHERE ps.id=$id");
		
		$sample_name = $sample_name_parts[0]["name"]."_0".$sample_name_parts[0]["process_id"];
		$folder_override = $base_folder_path."Sample_".$sample_name."/";
		$db->executeStmt("UPDATE processed_sample SET folder_override='$folder_override' WHERE id=$id");
	}
}

start_test($name);

//init
$db = DB::getInstance("NGSD_TEST");

//merge fastq germline sample:
reset_test_db($db, $name);
clear_analysis_folder();

create_folder("Sample_DNA220004_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220004_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220004_01_L001_R2_001.fastq.gz"], "dummy.txt");
create_folder("Sample_DNA220004_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"DNA220004_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220004_02_L001_R2_001.fastq.gz"], "dummy2.txt");

check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=4"), 1); // should be 1 before merging
check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220004_01 -into DNA220004_02");
check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=4"), 0); // should be 0 after merging

check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_01_L001_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_01_L001_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_02_L001_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_02_L001_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/dummy2.txt");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_DNA220004_01/");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_DNA220004_01/dummy.txt");

check($db->getValue("SELECT COUNT(*) FROM detected_variant WHERE processed_sample_id=4"), 0);
check($db->getValue("SELECT COUNT(*) FROM detected_variant WHERE processed_sample_id!=4"), 14);

check($db->getValue("SELECT COUNT(*) FROM cnv"), 24);
check($db->getValue("SELECT COUNT(*) FROM cnv_callset"), 2);

check($db->getValue("SELECT COUNT(*) FROM sv_deletion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_duplication"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_insertion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_inversion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_translocation"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_callset"), 2);

check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id=4"), 0);
check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id!=4"), 30);

check($db->getValue("SELECT COUNT(*) FROM merged_processed_samples WHERE processed_sample_id = 4 AND merged_into = 12"), 1); //merge correctly added

// merge fastq tumor sample with secondary analysis:

//reinit db
reset_test_db($db, $name);
//add secondary analysis here as gsvar path is dependend on the repository directory:
$db->executeStmt("INSERT INTO secondary_analysis (type, gsvar_file) VALUES ('somatic', '".data_folder()."/+analysis/merge_samples/Somatic_DNA220001_01-DNA220002_01/DNA220001_01-DNA220002_01.GSvar"."')");
clear_analysis_folder();

create_folder("Sample_DNA220001_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220001_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220001_01_L001_R2_001.fastq.gz"], "dummy.txt");
create_folder("Sample_DNA220001_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"DNA220001_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220001_02_L001_R2_001.fastq.gz"]);
create_folder("Somatic_DNA220001_01-DNA220002_01", [], "DNA220001_01-DNA220002_01.GSvar");

check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=1"), 1); // should be 1 before merging
check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220001_01 -into DNA220001_02");
check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=1"), 0); // should be 0 after merging

check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220001_02/DNA220001_01_L001_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220001_02/DNA220001_01_L001_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220001_02/DNA220001_02_L001_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220001_02/DNA220001_02_L001_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_DNA220001_01/");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_DNA220001_01/dummy.txt");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_analysis/Somatic_DNA220001_01-DNA220002_01/");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_analysis/Somatic_DNA220001_01-DNA220002_01/DNA220001_01-DNA220002_01.GSvar");

check($db->getValue("SELECT COUNT(*) FROM detected_somatic_variant WHERE processed_sample_id_tumor=1"), 0);
check($db->getValue("SELECT COUNT(*) FROM detected_somatic_variant WHERE processed_sample_id_tumor!=1"), 15);

check($db->getValue("SELECT COUNT(*) FROM somatic_cnv"), 26);
check($db->getValue("SELECT COUNT(*) FROM somatic_cnv_callset"), 2);

check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id=1"), 0);
check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id!=1"), 30);

check($db->getValue("SELECT COUNT(*) FROM secondary_analysis"), 0);

check($db->getValue("SELECT COUNT(*) FROM merged_processed_samples WHERE processed_sample_id = 1 AND merged_into = 9"), 1); //merge correctly added


//merge bam germline samples:
reset_test_db($db, $name);
clear_analysis_folder();

create_folder("Sample_DNA220004_01", [data_folder()."vc_strelka2_tu_in.bam"=>"DNA220004_01.bam"], "dummy.txt");
create_folder("Sample_DNA220004_02", [data_folder()."vc_strelka2_no_in.bam"=>"DNA220004_02.bam"]);

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220004_01 -into DNA220004_02 -sys ".data_folder()."merge_samples_system.txt");

check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_01_BamToFastq_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_01_BamToFastq_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_02_BamToFastq_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_02_BamToFastq_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_DNA220004_01/");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_DNA220004_01/dummy.txt");

check($db->getValue("SELECT COUNT(*) FROM detected_variant WHERE processed_sample_id=4"), 0);
check($db->getValue("SELECT COUNT(*) FROM detected_variant WHERE processed_sample_id!=4"), 14);

check($db->getValue("SELECT COUNT(*) FROM cnv"), 24);
check($db->getValue("SELECT COUNT(*) FROM cnv_callset"), 2);

check($db->getValue("SELECT COUNT(*) FROM sv_deletion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_duplication"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_insertion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_inversion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_translocation"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_callset"), 2);

check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id=4"), 0);
check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id!=4"), 30);

check($db->getValue("SELECT COUNT(*) FROM merged_processed_samples WHERE processed_sample_id = 4 AND merged_into = 12"), 1); //merge correctly added


//merge tumor RNA samples:
reset_test_db($db, $name);
clear_analysis_folder();

create_folder("Sample_RNA220003_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"RNA220003_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"RNA220003_01_L001_R2_001.fastq.gz", data_folder()."merge_samples_in1_L001_index_001.fastq.gz"=>"RNA220003_01_L001_index_001.fastq.gz"], "dummy.txt");
create_folder("Sample_RNA220003_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"RNA220003_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"RNA220003_02_L001_R2_001.fastq.gz", data_folder()."merge_samples_in2_L001_index_001.fastq.gz"=>"RNA220003_02_L001_index_001.fastq.gz"], "dummy2.txt");

check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=3"), 1); // should be 1 before merging
check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps RNA220003_01 -into RNA220003_02");
check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=3"), 0); // should be 0 after merging

check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_01_L001_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_01_L001_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_01_L001_index_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_02_L001_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_02_L001_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_02_L001_index_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/dummy2.txt");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_RNA220003_01/");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_RNA220003_01/dummy.txt");

check($db->getValue("SELECT COUNT(*) FROM expression WHERE processed_sample_id=3"), 0);
check($db->getValue("SELECT COUNT(*) FROM expression WHERE processed_sample_id!=3"), 4);

check($db->getValue("SELECT COUNT(*) FROM expression_exon WHERE processed_sample_id=3"), 0);
check($db->getValue("SELECT COUNT(*) FROM expression_exon WHERE processed_sample_id!=3"), 4);

check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id=3"), 0);
check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id!=3"), 30);

check($db->getValue("SELECT COUNT(*) FROM merged_processed_samples WHERE processed_sample_id = 3 AND merged_into = 11"), 1); //merge correctly added

// tumor sample: somatic report config exists -> has to fail
clear_analysis_folder();
create_folder("Sample_DNA220007_01", [data_folder()."vc_strelka2_tu_in.bam"=>"DNA220007_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220007_01_L001_R2_001.fastq.gz"]);
create_folder("Sample_DNA220007_02", [data_folder()."vc_strelka2_tu_in.bam"=>"DNA220007_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220007_02_L001_R2_001.fastq.gz"]);

check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=7"), 1); // should be 1 before merging
check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220007_01 -into DNA220007_02", false); // should fail because of somatic report config!
check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=7"), 1); // should still be 1 because merging failed

// normal sample: somatic report config exists -> is ok keep secondary analysis
clear_analysis_folder();
create_folder("Sample_DNA220008_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220008_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220008_01_L001_R2_001.fastq.gz"]);
create_folder("Sample_DNA220008_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"DNA220008_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220008_02_L001_R2_001.fastq.gz"]);

check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=8"), 1); // should be 1 before merging
check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220008_01 -into DNA220008_02");
check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=8"), 0); // should be 0 after merging


// germline sample: germline report config exists -> has to fail
clear_analysis_folder();
create_folder("Sample_DNA220006_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220006_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220006_01_L001_R2_001.fastq.gz"]);
create_folder("Sample_DNA220006_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"DNA220006_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220006_02_L001_R2_001.fastq.gz"]);

check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=6"), 1); // should be 1 before merging
check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220006_01 -into DNA220006_02", false); // should fail because of report config!
check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=6"), 1); // should still be 1 because merging failed

//merge cram germline samples: (existing index fastq files should be copied but cram files should still be turend into fastqs (index fastqs from demultiplexing with rna samples))
reset_test_db($db, $name);
clear_analysis_folder();

create_folder("Sample_DNA220004_01", [data_folder()."merge_samples.cram"=>"DNA220004_01.cram", data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220004_01_L001_index_001.fastq.gz"], "dummy.txt");
create_folder("Sample_DNA220004_02", [data_folder()."merge_samples.cram"=>"DNA220004_02.cram", data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220004_02_L001_index_001.fastq.gz"]);

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220004_01 -into DNA220004_02 -sys ".data_folder()."merge_samples_system.txt");
check($db->getValue("SELECT scheduled_for_resequencing FROM processed_sample WHERE id=4"), 0); // should be 0 after

check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_01_BamToFastq_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_01_BamToFastq_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_01_L001_index_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_02_BamToFastq_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_02_BamToFastq_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_DNA220004_02/DNA220004_02_L001_index_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_DNA220004_01/");
check_file_exists(data_folder()."/+analysis/merge_samples/+merged_samples/Sample_DNA220004_01/dummy.txt");

check($db->getValue("SELECT COUNT(*) FROM detected_variant WHERE processed_sample_id=4"), 0);
check($db->getValue("SELECT COUNT(*) FROM detected_variant WHERE processed_sample_id!=4"), 14);

check($db->getValue("SELECT COUNT(*) FROM cnv"), 24);
check($db->getValue("SELECT COUNT(*) FROM cnv_callset"), 2);

check($db->getValue("SELECT COUNT(*) FROM sv_deletion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_duplication"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_insertion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_inversion"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_translocation"), 5);
check($db->getValue("SELECT COUNT(*) FROM sv_callset"), 2);

check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id=4"), 0);
check($db->getValue("SELECT COUNT(*) FROM processed_sample_qc WHERE processed_sample_id!=4"), 30);

check($db->getValue("SELECT COUNT(*) FROM merged_processed_samples WHERE processed_sample_id = 4 AND merged_into = 12"), 1); //merge correctly added


//clean up rest up copied data
clear_analysis_folder();
end_test();


?>