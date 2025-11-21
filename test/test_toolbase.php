<?php

include("framework.php");

//##################################################################################
start_test("ToolBase::storeTDX");

//create dummy tool
$tool = new ToolBase("test_toolbase", "0.0a", "full description");

//extract version
$ver = $tool->extractVersion("cut");
check($ver!="n/a", true);

//exec
list($stdout, $stderr, $exit_code) = $tool->exec("echo", "bla && (echo bla2 >&2) && exit 17", false, false, false);
check(count($stdout)==1, true);
check(implode("", $stdout)=="bla", true);
check(count($stderr)==1, true);
check(implode("", $stderr)=="bla2", true);
check($exit_code==17, true);

end_test();

//exec parallel
start_test("ToolBase::execParallel");
//Test1
$log_file1 = output_folder()."/execParallel_1.log";
$tool->setLogFile($log_file1);
$out_file1 = output_folder()."/execParallel_command1.out";
$out_file2 = output_folder()."/execParallel_command2.out";
$out_file3 = output_folder()."/execParallel_command3.out";
$out_file4 = output_folder()."/execParallel_command4.out";
$out_file5 = output_folder()."/execParallel_command5.out";

$jobs = [
	["command1_exit0", "sleep 10s; echo 'command 1 finished!' | tee {$out_file1}"],
	["command2_exit0", "sleep 10s; echo 'command 2 finished!' | tee {$out_file2}"],
	["command3_exit0", "sleep 10s; echo 'command 3 finished!' | tee {$out_file3}"],
	["command4_exit0", "sleep 10s; echo 'command 4 finished!' | tee {$out_file4}"],
	["command5_exit0", "sleep 10s; echo 'command 5 finished!' | tee {$out_file5}"]
];

$tool->execParallel($jobs, 10);

check(file_exists($out_file1), true);
check(file($out_file1), ["command 1 finished!\n"]);

check(file_exists($out_file2), true);
check(file($out_file2), ["command 2 finished!\n"]);

check(file_exists($out_file3), true);
check(file($out_file3), ["command 3 finished!\n"]);

check(file_exists($out_file4), true);
check(file($out_file4), ["command 4 finished!\n"]);

check(file_exists($out_file5), true);
check(file($out_file5), ["command 5 finished!\n"]);

//Test2

$log_file2 = output_folder()."/execParallel_2.log";
$tool->setLogFile($log_file2);
$out_file6 = output_folder()."/execParallel_command6.out";
$out_file7 = output_folder()."/execParallel_command7.out";
$out_file8 = output_folder()."/execParallel_command8.out";

$jobs = [
	["command6_exit0", "sleep 10s; echo 'command 6 finished!' | tee {$out_file6}"],
	["command6_exit123", "sleep 10s; exit 123"],
	["command7_exit0", "sleep 10s; echo 'command 7 finished!' | tee {$out_file7}"],
	["command7_exit0_stdout", "sleep 10s; echo 'Totally random string to be logged in Tool output logfile!'"],
	["command8_exit0", "sleep 10s; echo 'command 8 finished!' | tee {$out_file8}"]
];

$tool->execParallel($jobs, 10, true, false, true);

check(file_exists($out_file6), true);
check(file($out_file6), ["command 6 finished!\n"]);

//job failed and gets logged in log file:
$log_file2_buffer = implode("\n", file($log_file2));
check(str_contains($log_file2_buffer, "WARNING: 'Processing of job command6_exit123_"), true);

check(file_exists($out_file7), true);
check(file($out_file7), ["command 7 finished!\n"]);

//job only reports to STDOUT:
check(str_contains($log_file2_buffer, "Totally random string to be logged in Tool output logfile!"), true);

check(file_exists($out_file8), true);
check(file($out_file8), ["command 8 finished!\n"]);



end_test();

?>
