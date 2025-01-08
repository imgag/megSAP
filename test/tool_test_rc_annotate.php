<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_annotate";
$gtf_file = get_path("data_folder")."/dbs/gene_annotations/GRCh38.gtf";

start_test($name);

// Annotate Ensembl gene identifiers with HGNC gene identifier/gene names
$in_file1 = data_folder()."rc_annotate_in1.tsv";
$out_file1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/Tools/{$name}.php -in {$in_file1} -out {$out_file1} -annotationIds gene_name,gene_biotype -gtf_file {$gtf_file}");
check_file($out_file1, data_folder().$name."_out1.tsv");

end_test();

?>