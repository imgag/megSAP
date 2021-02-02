<?php

require_once("framework.php");

$name = "create_circos_plot";
start_test($name);


//test 1: ClinCNV
// copy data 
$baf_file = output_folder()."clincnvTest_bafs.igv";
$seg_file = output_folder()."clincnvTest_cnvs_clincnv.seg";
$cnv_file = output_folder()."clincnvTest_cnvs_clincnv.tsv";
$roh_file = output_folder()."clincnvTest_rohs.tsv";
$sv_file = output_folder()."clincnvTest_manta_var_structural.bedpe";
copy(data_folder()."create_circos_plot_in_baf.igv", $baf_file);
copy(data_folder()."create_circos_plot_in_seg_clincnv.seg", $seg_file);
copy(data_folder()."create_circos_plot_in_cnv_clincnv.tsv", $cnv_file);
copy(data_folder()."create_circos_plot_in_rohs.tsv", $roh_file);
copy(data_folder()."create_circos_plot_in_sv.bedpe", $sv_file);

// create plot
check_exec("php ".src_folder()."/NGS/".$name.".php -folder ".output_folder()." -name clincnvTest");
check_file(output_folder()."clincnvTest_circos.png", data_folder().$name."_out_circos_clincnv.png");

end_test();

?>