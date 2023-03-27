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

start_test($name);

//init
$db = DB::getInstance("NGSD_TEST");

//merge fastq germline sample:
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
clear_analysis_folder();

create_folder("Sample_DNA220004_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220004_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220004_01_L001_R2_001.fastq.gz"], "dummy.txt");
create_folder("Sample_DNA220004_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"DNA220004_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220004_02_L001_R2_001.fastq.gz"], "dummy2.txt");

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220004_01 -into DNA220004_02");

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
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
//add secondary analysis here as gsvar path is dependend on the repository directory:
$db->executeStmt("INSERT INTO secondary_analysis (type, gsvar_file) VALUES ('somatic', '".data_folder()."/+analysis/merge_samples/Somatic_DNA220001_01-DNA220002_01/DNA220001_01-DNA220002_01.GSvar"."')");
clear_analysis_folder();

create_folder("Sample_DNA220001_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220001_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220001_01_L001_R2_001.fastq.gz"], "dummy.txt");
create_folder("Sample_DNA220001_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"DNA220001_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220001_02_L001_R2_001.fastq.gz"]);
create_folder("Somatic_DNA220001_01-DNA220002_01", [], "DNA220001_01-DNA220002_01.GSvar");

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220001_01 -into DNA220001_02");

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
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
clear_analysis_folder();

create_folder("Sample_DNA220004_01", [data_folder()."vc_strelka2_tu_in.bam"=>"DNA220004_01.bam"], "dummy.txt");
create_folder("Sample_DNA220004_02", [data_folder()."vc_strelka2_no_in.bam"=>"DNA220004_02.bam"]);

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220004_01 -into DNA220004_02");

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
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
clear_analysis_folder();

create_folder("Sample_RNA220003_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"RNA220003_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"RNA220003_01_L001_R2_001.fastq.gz"], "dummy.txt");
create_folder("Sample_RNA220003_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"RNA220003_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"RNA220003_02_L001_R2_001.fastq.gz"], "dummy2.txt");

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps RNA220003_01 -into RNA220003_02");

check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_01_L001_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_01_L001_R2_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_02_L001_R1_001.fastq.gz");
check_file_exists(data_folder()."/+analysis/merge_samples/Sample_RNA220003_02/RNA220003_02_L001_R2_001.fastq.gz");
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

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220007_01 -into DNA220007_02", false); // should fail because of somatic report config!

// normal sample: somatic report config exists -> is ok keep secondary analysis
clear_analysis_folder();
create_folder("Sample_DNA220008_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220008_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220008_01_L001_R2_001.fastq.gz"]);
create_folder("Sample_DNA220008_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"DNA220008_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220008_02_L001_R2_001.fastq.gz"]);

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220008_01 -into DNA220008_02");

// germline sample: germline report config exists -> has to fail
clear_analysis_folder();
create_folder("Sample_DNA220006_01", [data_folder()."merge_samples_in1_L001_R1_001.fastq.gz"=>"DNA220006_01_L001_R1_001.fastq.gz", data_folder()."merge_samples_in1_L001_R2_001.fastq.gz"=>"DNA220006_01_L001_R2_001.fastq.gz"]);
create_folder("Sample_DNA220006_02", [data_folder()."merge_samples_in2_L001_R1_001.fastq.gz"=>"DNA220006_02_L001_R1_001.fastq.gz", data_folder()."merge_samples_in2_L001_R2_001.fastq.gz"=>"DNA220006_02_L001_R2_001.fastq.gz"]);

check_exec("php ".src_folder()."/Tools/".$name.".php -db NGSD_TEST -ps DNA220006_01 -into DNA220006_02", false); // should fail because of report config!

//clean up rest up copied data
clear_analysis_folder();
// */
end_test();


?>