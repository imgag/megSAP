<?php

require_once("framework.php");

$name = "an_vep";
start_test($name);

//standard
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in1.vcf -out $out_file1 --log ".output_folder().$name."_out1.log");
remove_lines_containing($out_file1, array("##VEP=\"v"));
check_file($out_file1, data_folder().$name."_out1.vcf", true);

//empty input VCF
$out_file2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -out $out_file2 --log ".output_folder().$name."_out2.log");
remove_lines_containing($out_file2, array("##VEP=\"v"));
check_file($out_file2, data_folder().$name."_out2.vcf", true);

//somatic
$out_file3 = output_folder().$name."_out3.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -somatic -in ".data_folder().$name."_in1.vcf -out $out_file3 --log ".output_folder().$name."_out3.log");
remove_lines_containing($out_file3, array("##VEP=\"v"));
check_file($out_file3, data_folder().$name."_out3.vcf", true);

//NA12878_38 head (disease group test)
$out_file4 = output_folder().$name."_out4.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in3.vcf -out $out_file4 --log ".output_folder().$name."_out4.log");
remove_lines_containing($out_file4, array("##VEP=\"v"));
check_file($out_file4, data_folder().$name."_out4.vcf", true);

//NA12878_38 head zipped (disease group test)
$out_file5 = output_folder().$name."_out5.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in4.vcf.gz -out $out_file5 --log ".output_folder().$name."_out5.log");
remove_lines_containing($out_file5, array("##VEP=\"v"));
check_file($out_file5, data_folder().$name."_out5.vcf", true);

end_test();

?>
