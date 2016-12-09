<?php

require_once("framework.php");

$name = "mip_generator";

start_test($name);

$in_file = data_folder()."/".$name."_in1.bed";

$out_file1 = output_folder().$name."_out1_mips_picked.txt";
check_exec("php ".src_folder()."/Tools/".$name.".php -target $in_file -project_name ${name}_out1 -mode blood -out_folder ".output_folder());
check_file($out_file1, data_folder().$name."_out1.picked_mips.txt");

$out_file2 = output_folder().$name."_out2_mips_picked.txt";
check_exec("php ".src_folder()."/Tools/".$name.".php -target $in_file -project_name ${name}_out2 -mode ffpe -out_folder ".output_folder());
check_file($out_file2, data_folder().$name."_out2.picked_mips.txt");

end_test();

?>