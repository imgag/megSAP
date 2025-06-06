<?php

require_once("framework.php");

$name = "vc_manta";
start_test($name);

########################## germline (exome) ##################
$out_file1 = output_folder().$name."_out1.vcf.gz";
$manta_evidence_dir = output_folder()."/manta_evid_1";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in1.bam -out $out_file1 -evid_dir $manta_evidence_dir -exome -target ".data_folder().$name."_enrichment.bed -threads 2 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf.gz");

//TODO: change files to GRCh38
########################## tumor-only ########################
$out_file2 = output_folder().$name."_out2.vcf.gz";
$manta_evidence_dir = output_folder()."/manta_evid_2";
check_exec("php ".src_folder()."/Tools/{$name}.php -t_bam ".data_folder()."vc_strelka2_tu_in.bam -out $out_file2 -evid_dir $manta_evidence_dir -threads 2 --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.vcf.gz");

//TODO: change files to GRCh38
########################## somatic ###########################
$out_file3 = output_folder().$name."_out3.vcf.gz";
$manta_evidence_dir = output_folder()."/manta_evid_3";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder()."vc_strelka2_no_in.bam -t_bam ".data_folder()."vc_strelka2_tu_in.bam -out $out_file3 -evid_dir $manta_evidence_dir -threads 2 --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.vcf.gz");

########################## RNA ###############################
$out_file4 = output_folder().$name."_out4.vcf.gz";
$manta_evidence_dir = output_folder()."/manta_evid_4";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in4.bam -out $out_file4 -evid_dir $manta_evidence_dir -rna -threads 2 --log ".output_folder().$name."_out4.log");
check_file($out_file4, data_folder().$name."_out4.vcf.gz");

########################## INV fix ###############################
$out_file5 = output_folder().$name."_out5.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in5.bam -out $out_file5 -threads 2 --log ".output_folder().$name."_out5.log");
check_file($out_file5, data_folder().$name."_out5.vcf.gz");
$out_file6 = output_folder().$name."_out6.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in5.bam -skip_inv_merging -out $out_file6 -threads 2 --log ".output_folder().$name."_out6.log");
check_file($out_file6, data_folder().$name."_out6.vcf.gz");

end_test();

?>


