<?php

require_once("framework.php");

$name = "create_methyl_plot";
start_test($name);

//create cohort folder
exec2("mkdir -p  ".output_folder()."/cohort/ProcessingSystemForTesting");

//test 1 bam
exec2("cp ".data_folder().$name."_in1.bam ".output_folder()."Testsample_01.bam");
exec2("cp ".data_folder().$name."_in1.bam.bai ".output_folder()."Testsample_01.bam.bai");
exec2("cp ".data_folder().$name."_in1_phasing_track.bed ".output_folder()."Testsample_01_phasing_track.bed");

check_exec("php ".src_folder()."/Tools/".$name.".php -folder ".output_folder()." -test -name Testsample_01 -out ".output_folder().$name."_out1.tsv -custom_cohort_table ".data_folder().$name."_custom_cohort.tsv -regions ".data_folder().$name."_regions.tsv --log ".output_folder().$name."_1.log");
check_file(output_folder().$name."_out1.tsv", data_folder().$name."_out1.tsv");

check(filesize(output_folder()."methylartist/Testsample_01_GRB10_alt-TSS-DMR.png"), 321023, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_PLAGL1_alt-TSS-DMR.png"), 315217, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_PEG10_TSS-DMR.png"), 272116, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_PHF6.png"), 222518, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_RP2.png"), 264019, 2048);

check(filesize(output_folder()."methylartist/Testsample_01_GRB10_alt-TSS-DMR_cohort.png"), 180102, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_PLAGL1_alt-TSS-DMR_cohort.png"), 114136, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_PEG10_TSS-DMR_cohort.png"), 65202, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_PHF6_cohort.png"), 45832, 2048);
check(filesize(output_folder()."methylartist/Testsample_01_RP2_cohort.png"), 44453, 2048);


//test 2 cram
exec2("cp ".data_folder().$name."_in1.cram ".output_folder()."Testsample_02.cram");
exec2("cp ".data_folder().$name."_in1.cram.crai ".output_folder()."Testsample_02.cram.crai");
exec2("cp ".data_folder().$name."_in1_phasing_track.bed ".output_folder()."Testsample_02_phasing_track.bed");

check_exec("php ".src_folder()."/Tools/".$name.".php -folder ".output_folder()." -test -name Testsample_02 -out ".output_folder().$name."_out2.tsv -export_methylation_data -custom_cohort_table ".data_folder().$name."_custom_cohort.tsv -regions ".data_folder().$name."_regions.tsv -skip_plot --log ".output_folder().$name."_2.log");
check_file(output_folder().$name."_out2.tsv", data_folder().$name."_out2.tsv");
check_file(output_folder().$name."_out2_raw_data.tsv", data_folder().$name."_out2_raw_data.tsv");

check(file_exists(output_folder()."methylartist/Testsample_02_GRB10_alt-TSS-DMR.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_PLAGL1_alt-TSS-DMR.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_PEG10_TSS-DMR.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_PHF6.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_RP2.png"), false);

check(file_exists(output_folder()."methylartist/Testsample_02_GRB10_alt-TSS-DMR_cohort.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_PLAGL1_alt-TSS-DMR_cohort.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_PEG10_TSS-DMR_cohort.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_PHF6_cohort.png"), false);
check(file_exists(output_folder()."methylartist/Testsample_02_RP2_cohort.png"), false);

end_test();

?>