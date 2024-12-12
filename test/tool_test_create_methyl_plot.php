<?php

require_once("framework.php");

$name = "create_methyl_plot";
start_test($name);



//test 1 bam
exec2("cp ".data_folder().$name."_in1.bam ".output_folder()."Testsample_01.bam");
exec2("cp ".data_folder().$name."_in1.bam.bai ".output_folder()."Testsample_01.bam.bai");

check_exec("php ".src_folder()."/NGS/".$name.".php -folder ".output_folder()." -name Testsample_01 -out ".output_folder().$name."_out1.tsv -regions ".data_folder().$name."_regions.tsv --log ".output_folder().$name."_1.log");
check_file(output_folder().$name."_out1.tsv", data_folder().$name."_out1.tsv");

check(filesize(output_folder()."methylartist/Testsample_01_GRB10_alt-TSS-DMR.png"), 255008, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_PLAGL1_alt-TSS-DMR.png"), 255559, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_PEG10_TSS-DMR.png"), 258754, 2048);


//test 2 cram
exec2("cp ".data_folder().$name."_in1.cram ".output_folder()."Testsample_02.cram");
exec2("cp ".data_folder().$name."_in1.cram.crai ".output_folder()."Testsample_02.cram.crai");

check_exec("php ".src_folder()."/NGS/".$name.".php -folder ".output_folder()." -name Testsample_02 -out ".output_folder().$name."_out2.tsv -regions ".data_folder().$name."_regions.tsv -skip_plot --log ".output_folder().$name."_2.log");
check_file(output_folder().$name."_out2.tsv", data_folder().$name."_out1.tsv");

check(file_exists(output_folder()."methylartist/Testsample_02_GRB10_alt-TSS-DMR.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_PLAGL1_alt-TSS-DMR.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_PEG10_TSS-DMR.png"), false);


end_test();

?>