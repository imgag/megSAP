<?php

require_once("framework.php");

$name = "an_vep";
start_test($name);

//standard
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in1.vcf -out $out_file1 --log ".output_folder().$name."_out1.log");
remove_lines_containing($out_file1, array("##VEP=\"v", "##VEP-command-line='", "##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file1, data_folder().$name."_out1.vcf", true);

//gzipped input VCF
$out_file2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in2.vcf.gz -out $out_file2 --log ".output_folder().$name."_out2.log");
remove_lines_containing($out_file2, array("##VEP=\"v", "##VEP-command-line='", "##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file2, data_folder().$name."_out2.vcf", true);

//empty input VCF
$out_file_empty = output_folder().$name."_out_empty.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -out $out_file_empty --log ".output_folder().$name."_out_empty.log");
remove_lines_containing($out_file_empty, array("##VEP=\"v", "##VEP-command-line='", "##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file_empty, data_folder().$name."_out_empty.vcf", true);

//variants to score with splicing variants
$out_file3 = output_folder().$name."_out3.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in3.vcf -out $out_file3 --log ".output_folder().$name."_out3.log");
remove_lines_containing($out_file3, array("##VEP=\"v", "##VEP-command-line='", "##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file3, data_folder().$name."_out3.vcf", true);

//DRAGEN
$out_file_dragen = output_folder().$name."_out_dragen.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in_dragen.vcf -out $out_file_dragen --log ".output_folder().$name."_out_dragen.log");
remove_lines_containing($out_file_dragen, array("##VEP=\"v", "##VEP-command-line='", "##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file_dragen, data_folder().$name."_out_dragen.vcf", true);

//tests with NGSD disease group annotations - these tests will be run only if NGSD is available (i.e. not in nightly tests)
if (db_is_enabled("NGSD")) 
{
	//standard
	$out_file_db1 = output_folder().$name."_out_db1.vcf";
	check_exec("php ".src_folder()."/NGS/{$name}.php -test -ps_name NA12878_03 -in ".data_folder().$name."_in1.vcf -out $out_file_db1 --log ".output_folder().$name."_out_db1.log");
	remove_lines_containing($out_file_db1, array("##VEP=\"v", "##VEP-command-line='", "##INFO=<ID=NGSD_GENE_INFO,"));
	check_file($out_file_db1, data_folder().$name."_out_db1.vcf", true);

	//somatic
	$out_file_som = output_folder().$name."_out_som.vcf";
	check_exec("php ".src_folder()."/NGS/{$name}.php -test -somatic -in ".data_folder().$name."_in1.vcf -out $out_file_som --log ".output_folder().$name."_out_som.log");
	remove_lines_containing($out_file_som, array("##VEP=\"v", "##VEP-command-line='", "##INFO=<ID=NGSD_GENE_INFO,"));
	check_file($out_file_som, data_folder().$name."_out_som.vcf", true);
}

//Somatic VICC data
$out_file_som_vicc = output_folder().$name."_out_som_vicc.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in5_som_vicc.vcf -out $out_file_som_vicc -somatic --log ".output_folder().$name."_out_som_vicc.log");
remove_lines_containing($out_file_som_vicc, array("##VEP=\"v", "##VEP-command-line='", "##INFO=<ID=NGSD_GENE_INFO,", "##INFO=<ID=NGSD_GROUP"));
check_file($out_file_som_vicc, data_folder().$name."_out_som_vicc.vcf", true);

end_test();

?>
