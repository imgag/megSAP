<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_annotate";
$build = "GRCh38";
$gtfFile = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

start_test($name);

// Annotate Ensembl gene identifiers with HGNC gene identifier/gene names
$in_file1 = data_folder()."rc_normalize_out1.tsv";
$out_file1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in {$in_file1} -out {$out_file1} -annotationIds gene_name -gtfFile {$gtfFile}");
check_file($out_file1, data_folder().$name."_out1.tsv");

end_test();

?>