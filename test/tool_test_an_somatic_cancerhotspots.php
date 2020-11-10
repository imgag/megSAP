<?php

require_once("framework.php");

$name = "an_somatic_cancerhotspots";

start_test($name);
check_exec("php ".src_folder()."/NGS/{$name}.php -in " . data_folder() ."{$name}_in1.vcf -out " . output_folder(). "{$name}_out1.vcf");
check_file(output_folder(). "{$name}_out1.vcf", data_folder().$name."_ref1.vcf");
end_test();

?>
