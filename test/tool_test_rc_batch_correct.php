<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/functions.php");

$name = "rc_batch_correct";

start_test($name);

$command_only = true;

// test normalization
$in_data = data_folder().$name."_in1.tsv.gz";
$in_batch = data_folder().$name."_batch1.tsv";
$out = output_folder().$name."_out1.tsv.gz";
$args = [];
$args[] = "-i {$in_data}";
$args[] = "-b {$in_batch}";
$args[] = "-o {$out}";
$apptainer_cmd = execApptainer("python", "python3", src_folder()."/Tools/{$name}.py ".implode(" ", $args), [src_folder()."/Tools/{$name}.py", $in_data, $in_batch], [output_folder()], $command_only);
check_exec($apptainer_cmd);
check_file($out, data_folder().$name."_out1.tsv.gz");


// test normalization with covariants
$in_data = data_folder().$name."_in2.tsv.gz";
$in_batch = data_folder().$name."_batch2.tsv";
$in_covar = data_folder().$name."_covar2.tsv";
$out = output_folder().$name."_out2.tsv.gz";
$args = [];
$args[] = "-i {$in_data}";
$args[] = "-b {$in_batch}";
$args[] = "-c {$in_covar}";
$args[] = "-o {$out}";
$apptainer_cmd = execApptainer("python", "python3", src_folder()."/Tools/{$name}.py ".implode(" ", $args), [src_folder()."/Tools/{$name}.py", $in_data, $in_batch, $in_covar], [dirname($out)], $command_only);
check_exec($apptainer_cmd);
check_file($out, data_folder().$name."_out2.tsv.gz");

end_test();

?>