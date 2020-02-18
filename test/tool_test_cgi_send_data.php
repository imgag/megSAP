<?php
require_once("framework.php");

$name = "cgi_send_data";

start_test($name);

//column names in CGI SNV file for checks
$col_names_mut = array("driver_statement","gene_role","transcript");

$out_zip = output_folder() ."/cgi_results1.zip";

$out_dir1 = output_folder() ."/out1/";
$out_dir2 = output_folder() ."/out2/";
$out_dir3 = output_folder() ."/out3/";
$out_dir4 = output_folder() ."/out4/";
exec2("mkdir $out_dir1 $out_dir2 $out_dir3 $out_dir4");


//Case 1: input is vcf.gz file
$in_file_snv1 = data_folder() . "/cgisenddata_in1_var_annotated.vcf.gz";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv1 -cancertype SK -out $out_zip");
exec2("unzip -n $out_zip -d " . $out_dir1 );
//check mutation file for columns which are later annotated to GSvar (sometimes CGI adds/removes columns, order changes often)
check_column_exists($out_dir1."mutation_analysis.tsv",$col_names_mut);
check_file_exists($out_dir1."drug_prescription.tsv");
check_file_exists($out_dir1."drug_prescription_bioactivities.tsv");


//Case 2: input file is .vcf file (output must be the same as if vcf.gz as input)
$in_file_snv2 = data_folder() . "/cgisenddata_in2_var_annotated.vcf";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv2 -cancertype SK -out $out_zip");
exec2("unzip -n $out_zip -d " . $out_dir2 );
check_column_exists($out_dir2."mutation_analysis.tsv",$col_names_mut);


//Case 3: input file .vcf and ClinCNV CNV
$in_file_snv3 = data_folder() . "/cgisenddata_in3_var_annotated.vcf";
$in_file_cnv3 = data_folder() . "/cgisenddata_in3_clincnv.tsv";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv3 -cnas $in_file_cnv3 -cancertype SK -out $out_zip");
exec2("unzip -n $out_zip -d " . $out_dir3 );
check_column_exists($out_dir3."mutation_analysis.tsv",$col_names_mut);
check_tsv_file($out_dir3."cna_analysis.tsv",data_folder()."cgisenddata_ref3_cgi_cnv_analysis.tsv");
check_file_exists($out_dir3."drug_prescription.tsv");
check_file_exists($out_dir3."drug_prescription_bioactivities.tsv");


//Case 4: Upload CNV file with target region set
$in_file_cnv4 = data_folder() . "/cgisenddata_in3_clincnv.tsv";
$target_region = data_folder() . "/cgisenddata_target_region_genes.txt";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -cnas $in_file_cnv4 -cancertype SK -out $out_zip -t_region $target_region",true);
exec2("unzip -n $out_zip -d " . $out_dir4 );
check_tsv_file($out_dir4."cna_analysis.tsv",data_folder()."cgisenddata_ref4_cgi_cnv_analysis.tsv");


end_test();

?>