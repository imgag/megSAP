<?php

require_once("framework.php");

$name = "gsvar2vcf";
start_test($name);

$in_file = data_folder().$name."_in.tsv";
$out_file = output_folder().$name."_out.vcf";

check_exec("php ".src_folder()."/Auxilary/gsvar2vcf.php -in '$in_file' -out '$out_file'");
remove_lines_containing($out_file, array("##fileDate="));
check_file($out_file, data_folder().$name."_out.vcf");

end_test();

?>