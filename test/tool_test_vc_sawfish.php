<?php

require_once("framework.php");

$name = "vc_sawfish";
start_test($name);

//NOTE: test data generated from (same file as PacBio DeeepVariant tool test)

########################## Single ##################
$out_file1 = output_folder().$name."_out1.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -bams ".data_folder()."{$name}_in.bam -keep_common_cnvs -small_var_vcfs ".data_folder()."{$name}_in.vcf.gz -genders 'male' -out {$out_file1} -test -threads 4 --log ".output_folder()."{$name}_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf.gz");


#TODO
########################## Trio ##################

end_test();

?>


