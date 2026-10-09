<?php

require_once("framework.php");

//NOTE: Test data was generated from mito region from benchmark samples NA24385-lrGS-DEFAULT_01 (10% downsample), NA24385-lrGS-PacBio_01

$name = "vc_himito";
start_test($name);

//ONT
$out_file1 = output_folder().$name."_out1.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -data_type ont-r10 -threads 4 -name NA24385-ONT_01 -bam ".data_folder().$name."_in1.bam -out $out_file1  --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf.gz");

//PacBio
$out_file2 = output_folder().$name."_out2.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -data_type pacbio -threads 4 -name NA24385-PacBio_01 -bam ".data_folder().$name."_in2.cram -out $out_file2  --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.vcf.gz");

end_test();

?>


