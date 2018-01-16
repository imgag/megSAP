<?php
require_once("framework.php");

$name = "cgi_send_data";

start_test($name);

$out_folder1 = output_folder();
$in_file_snv1 = data_folder() . "/cgisenddata_in1_var_annotated.vcf.gz";
$in_file_cnv1 = data_folder() . "/cgisenddata_in1_cnvs.tsv";
$cancertype = "SK";

//column names in CGI SNV file for checks
$col_names_mut = array("driver_statement","gene_role","transcript");

$out_file_snv_analysis_1 = output_folder()."/cgisenddata_in1_cgi_mutation_analysis.tsv";


//standard: input is vcf.gz file
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv1 -cnas $in_file_cnv1 -cancertype $cancertype -o_folder $out_folder1");
check_tsv_file(output_folder()."cgisenddata_in1_cgi_cnv_analysis.tsv",data_folder()."cgisenddata_ref1_cgi_cnv_analysis.tsv");
//check mutation file for columns which are later annotated to GSvar (sometimes CGI adds/removes columns, order changes often)

check_column_exists($out_file_snv_analysis_1,$col_names_mut);
check_file_exists(output_folder()."cgisenddata_in1_cgi_drug_prescription.tsv");
check_file_exists(output_folder()."cgisenddata_in1_cgi_drug_prescription_bioactivities.tsv");


//input file is .vcf file, output must be the same as if vcf.gz as input
$in_file_snv2 = data_folder() . "/cgisenddata_in2_var_annotated.vcf";
$out_file_snv_analysis_2 = output_folder()."/cgisenddata_in2_cgi_mutation_analysis.tsv";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv2 -cancertype SK -o_folder $out_folder1");
check_column_exists($out_file_snv_analysis_2,$col_names_mut);


//input file .vcf and .cnv, where CNV file is empty
$in_file_snv3 = data_folder() . "/cgisenddata_in3_var_annotated.vcf";
$in_file_cnv3 = data_folder() . "/cgisenddata_in3_cnvs.tsv";
$out_file_snv_analysis_3 = output_folder()."/cgisenddata_in3_cgi_mutation_analysis.tsv";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv3 -cnas $in_file_cnv3 -cancertype SK -o_folder $out_folder1");
check_column_exists($out_file_snv_analysis_3,$col_names_mut);
check_file_exists(output_folder()."cgisenddata_in3_cgi_drug_prescription.tsv");
check_file_exists(output_folder()."cgisenddata_in3_cgi_drug_prescription_bioactivities.tsv");


//input file .vcf empty, cnv contains variants
$in_file_snv4 = data_folder() . "/cgisenddata_in4_var_annotated.vcf";
$in_file_cnv4 = data_folder() . "/cgisenddata_in4_cnvs.tsv";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv4 -cnas $in_file_cnv4 -cancertype SK -o_folder $out_folder1",true);
check_tsv_file(output_folder()."cgisenddata_in4_cgi_cnv_analysis.tsv",data_folder()."cgisenddata_ref4_cgi_cnv_analysis.tsv");
check_file_exists(output_folder()."cgisenddata_in4_cgi_drug_prescription.tsv");
check_file_exists(output_folder()."cgisenddata_in4_cgi_drug_prescription_bioactivities.tsv");


//input file cnv:  in_file_cnv4, snvs empty and target_region is set
$target_region = data_folder() . "/cgisenddata_target_region_genes.txt";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -cnas $in_file_cnv4 -cancertype SK -o_folder $out_folder1 -t_region $target_region",true);
check_tsv_file(output_folder()."cgisenddata_in4_cgi_cnv_analysis.tsv",data_folder()."cgisenddata_ref5_cgi_cnv_analysis.tsv");

?>