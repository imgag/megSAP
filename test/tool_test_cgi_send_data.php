<?php
require_once("framework.php");

$name = "cgi_send_data";

start_test($name);

//column names in CGI SNV file for checks
$col_names_mut = array("driver_statement","gene_role","transcript");

$out_zip = output_folder() ."/cgi_results1.zip";

$out_dir1 = output_folder() ."/out1/";
$out_dir2 = output_folder() ."/out2/";
exec2("mkdir $out_dir1 $out_dir2");

//Case 1: upload SNVS and small INDELS
$in_file_snv1 = data_folder() . "/cgisenddata_in1.GSvar";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -mutations $in_file_snv1 -cancertype SK -out $out_zip");
exec2("unzip -n $out_zip -d " . $out_dir1 );
//check mutation file for columns which are later annotated to GSvar (sometimes CGI adds/removes columns, order changes often)
check_column_exists($out_dir1."mutation_analysis.tsv",$col_names_mut);
check_file_exists($out_dir1."drug_prescription.tsv");
check_file_exists($out_dir1."drug_prescription_bioactivities.tsv");

//Case 2: Upload CNV file with target region set
$in_file_cnv2 = data_folder() . "/cgisenddata_in2_clincnv.tsv";
$target_region = data_folder() . "/cgisenddata_target_region_genes.txt";
check_exec("php ".src_folder()."/NGS/cgi_send_data.php -cnas $in_file_cnv2 -cancertype SK -out $out_zip -t_region $target_region",true);
exec2("unzip -n $out_zip -d " . $out_dir2 );
check_tsv_file($out_dir2."cna_analysis.tsv",data_folder()."cgisenddata_ref2_cgi_cnv_analysis.tsv");


end_test();

?>