<?php

require_once("framework.php");

$name = "backup_queue";
$file = data_folder().$name;

start_test($name);

//init NGSD
init_ngsd($name);

$input_folder = output_folder()."/input_folder/";
mkdir($input_folder);

//test single
$log_file = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/IMGAG/{$name}.php -test -in $input_folder -mode run -email ahtesto1 --log $log_file");
exec("grep \"Command: \" $log_file", $command_line);
check(strpos($command_line[0], " test@med.uni-tuebingen.de "), true); //email correct from DB
check(strpos($command_line[0], "/backup_run.php "), true); //mode correct
$command_line = [];

$log_file = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/IMGAG/{$name}.php -test -in $input_folder -mode run -include_raw_signal -email test@test.de --log $log_file");
exec("grep \"Command: \" $log_file", $command_line);
check(strpos($command_line[0], " test@test.de "), true); //email correct
check(strpos($command_line[0], "/backup_run.php "), true); //mode correct
check(strpos($command_line[0], " -include_raw_signal "), true); //argument correct
$command_line = [];


$log_file = output_folder().$name."_out3.log";
check_exec("php ".src_folder()."/IMGAG/{$name}.php -test -in $input_folder -mode user -email test@test.de --log $log_file");
exec("grep \"Command: \" $log_file", $command_line);
check(strpos($command_line[0], "/backup_user.php "), true); //mode correct
$command_line = [];


$log_file = output_folder().$name."_out4.log";
check_exec("php ".src_folder()."/IMGAG/{$name}.php -test -in $input_folder -mode project -email test@test.de --log $log_file");
exec("grep \"Command: \" $log_file", $command_line);
check(strpos($command_line[0], "/backup_project.php "), true); //mode correct
$command_line = [];

end_test();

?>
