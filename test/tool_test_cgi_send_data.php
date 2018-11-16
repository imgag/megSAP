<?php
require_once("framework.php");

$name = "cgi_send_data";

start_test($name);

//column names in CGI SNV file for checks
$col_names_mut = array("driver_statement","gene_role","transcript");

$out_zip = output_folder() ."/cgi_results1.zip";

//Case 1: input is vcf.gz file, CNV is CNVHunter file
$in_file_snv1 = data_folder() . "/cgisenddata_in1_var_annotated.vcf.gz";
$in_file_cnv1 = data_folder() . "/cgisenddata_in1_cnvs.tsv";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv1 -cnas $in_file_cnv1 -cancertype SK -out $out_zip");
exec2("unzip -n $out_zip -d " . output_folder() );
check_tsv_file(output_folder()."cna_analysis.tsv",data_folder()."cgisenddata_ref1_cgi_cnv_analysis.tsv");
//check mutation file for columns which are later annotated to GSvar (sometimes CGI adds/removes columns, order changes often)
check_column_exists(output_folder()."mutation_analysis.tsv",$col_names_mut);
check_file_exists(output_folder()."drug_prescription.tsv");
check_file_exists(output_folder()."drug_prescription_bioactivities.tsv");
clear_output_folder();

//Case 2: input file is .vcf file (output must be the same as if vcf.gz as input)
$in_file_snv2 = data_folder() . "/cgisenddata_in2_var_annotated.vcf";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv2 -cancertype SK -out $out_zip");
exec2("unzip -n $out_zip -d " . output_folder() );
check_column_exists(output_folder()."mutation_analysis.tsv",$col_names_mut);
clear_output_folder();

//Case 3: input file .vcf and ClinCNV CNV
$in_file_snv3 = data_folder() . "/cgisenddata_in3_var_annotated.vcf";
$in_file_cnv3 = data_folder() . "/cgisenddata_in3_clincnv.tsv";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv3 -cnas $in_file_cnv3 -cancertype SK -out $out_zip");
exec2("unzip -n $out_zip -d " . output_folder() );
check_column_exists(output_folder()."mutation_analysis.tsv",$col_names_mut);
check_tsv_file(output_folder()."cna_analysis.tsv",data_folder()."cgisenddata_ref3_cgi_cnv_analysis.tsv");
check_file_exists(output_folder()."drug_prescription.tsv");
check_file_exists(output_folder()."drug_prescription_bioactivities.tsv");
clear_output_folder();

//Case 4: Upload CNV file with target region set
$in_file_cnv4 = data_folder() . "/cgisenddata_in3_clincnv.tsv";
$target_region = data_folder() . "/cgisenddata_target_region_genes.txt";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -cnas $in_file_cnv4 -cancertype SK -out $out_zip -t_region $target_region",true);
exec2("unzip -n $out_zip -d " . output_folder() );
check_tsv_file(output_folder()."cna_analysis.tsv",data_folder()."cgisenddata_ref4_cgi_cnv_analysis.tsv");
clear_output_folder();

end_test();

?>