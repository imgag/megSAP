<?php

require_once("framework.php");

$name = "vc_ballele";
start_test($name);

$t_bam = data_folder()."vc_strelka2_tu_in.bam";
$n_bam = data_folder()."vc_strelka2_no_in.bam";

$snp_db = data_folder()."/".$name."_snps.vcf.gz";

// test1: somatic, with depth
// expected output: variants in target regions, two samples/columns
$out1 = output_folder().$name."_out1.igv";
check_exec("php ".src_folder()."/NGS/vc_ballele.php -bam $t_bam -n_bam $n_bam -out $out1 -vcf $snp_db -depth --log ".output_folder().$name."_1.log");
check_file($out1, data_folder().$name."_out1.igv");

// test2: germline
// expected output: variants in target regions, one sample/column
$out2 = output_folder().$name."_out2.igv";
check_exec("php ".src_folder()."/NGS/vc_ballele.php -bam $t_bam -out $out2 -vcf $snp_db --log ".output_folder().$name."_2.log");
check_file($out2, data_folder().$name."_out2.igv");

// test3: germline, with downsampling => WGS
// expected output: 3 variants
$out3 = output_folder().$name."_out3.igv";
check_exec("php ".src_folder()."/NGS/vc_ballele.php -bam $t_bam -out $out3 -vcf $snp_db -downsample -factor 10 --log ".output_folder().$name."_3.log");
check_file($out3, data_folder().$name."_out3.igv");


end_test();

?>