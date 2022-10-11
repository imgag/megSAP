<?php

require_once("framework.php");

$name = "baf_germline";
start_test($name);


// test1: strelka input
$vcf = data_folder()."/{$name}_in1.vcf.gz";
$out1 = output_folder()."{$name}_out1.igv";
check_exec("php ".src_folder()."/NGS/{$name}.php -name vc_strelka2_tu_in -out $out1 -vcf $vcf --log ".output_folder()."{$name}_1.log");
check_file($out1, data_folder()."{$name}_out1.igv");

// test3: dragen input
$vcf = data_folder()."/{$name}_in3.vcf.gz";
$out2 = output_folder()."{$name}_out2.igv";
check_exec("php ".src_folder()."/NGS/{$name}.php -name NA12878_03 -out $out2 -vcf $vcf --log ".output_folder()."{$name}_1.log");
check_file($out2, data_folder()."{$name}_out2.igv");

// test2: freebayes input, with downsampling
$vcf = data_folder()."/{$name}_in2.vcf.gz";
$out3 = output_folder()."{$name}_out3.igv";
check_exec("php ".src_folder()."/NGS/{$name}.php -name vc_strelka2_tu_in -out $out3 -vcf $vcf -downsample 10 --log ".output_folder()."{$name}_2.log");
check_file($out3, data_folder()."{$name}_out3.igv");

end_test();

?>