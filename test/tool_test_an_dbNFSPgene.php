<?php

require_once("framework.php");

$name = "an_dbNFSPgene";
start_test($name);

$out_file = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/NGS/".$name.".php -in ".data_folder().$name."_in1.GSvar -out $out_file");
check_file($out_file, data_folder().$name."_out1.GSvar");

end_test();

?>