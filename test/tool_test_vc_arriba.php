<?php

require_once("framework.php");

$name = "vc_arriba";
start_test($name);

$in_bam = data_folder().$name."_in1.bam";

$out_tsv = output_folder().$name."_out1.tsv";
$out_discarded = output_folder().$name."_out1.discarded.tsv";
$out_vcf = output_folder().$name."_out1.vcf";
$out_pdf = output_folder().$name."_out1.pdf";
$out_pic_dir = output_folder().$name."_out1";

check_exec("php ".src_folder()."/Tools/{$name}.php -bam $in_bam -out_fusions $out_tsv -out_discarded $out_discarded -out_vcf $out_vcf -out_pdf $out_pdf -out_pic_dir $out_pic_dir --log ".output_folder().$name."_out1.log");
#remove line as it contains the analysis date: '##AnalysisDate=%%%'
exec("sed -i '/##AnalysisDate/d' $out_tsv");
check_file($out_tsv, data_folder().$name."_ref1.tsv", true);

//diff vcf by hand as "check_file" expects valid vcf and arriba returns a "simplified vcf".
$logfile = $out_vcf."_diff";
exec("diff -b ".data_folder().$name."_ref1.vcf"." $out_vcf > $logfile 2>&1", $output, $return);
if($return != 0)
{
	print "Differenz in out vcfs. See $logfile \n";
}
check($return == 0, true);

end_test();

?>
