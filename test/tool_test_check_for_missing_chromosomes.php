<?php

require_once("framework.php");

$name = "check_for_missing_chromosomes";
$data_folder = data_folder();

start_test($name);

//VCF
$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in1.vcf");
check($output, []);

$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in2.vcf", false);
check_true(count($output)==1);
check_contains($output[0], "chr3, chr21");

//TSV
$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in3.tsv");
check($output, []);

$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in4.tsv", false);
check_true(count($output)==1);
check_contains($output[0], "chr3, chr10, chr19");

//GSvar
$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in5.GSvar");
check($output, []);

$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in6.GSvar", false);
check_true(count($output)==1);
check_contains($output[0], "chr10, chr11, chr12, chrX");

//VCF.GZ
$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in7.vcf.gz");
check($output, []);

$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in8.vcf.gz", false);
check_true(count($output)==1);
check_contains($output[0], "chr3, chr21");

//max_missing_perc
$output = check_exec("php ".src_folder()."/Tools/{$name}.php -in {$data_folder}/{$name}_in1.vcf -max_missing_perc 5", false);
check_true(count($output)==1);
check_contains($output[0], "chr1=92.6%");

end_test();


?>