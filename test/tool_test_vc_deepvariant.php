<?php

require_once("framework.php");

//NOTE: Test data was generated from NA12878_58: all reads overlapping the first 50 variants on chr22 (see make target 'vc_freebayes')

$name = "vc_deepvariant";
start_test($name);


########################## with target region ##########################

$out_file1 = output_folder().$name."_out1.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file1 -model_type WES -target ".data_folder().$name."_in.bed --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf.gz");

//with multiple threads - should give the same result as single-threaded
for ($i=2; $i<=4; ++$i)
{
	$out_file2 = output_folder().$name."_out2_{$i}threads.vcf.gz";
	check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file2 -model_type WES -target ".data_folder().$name."_in.bed -threads {$i} --log ".output_folder().$name."_out2_{$i}threads.log");
	check_file($out_file2, data_folder().$name."_out1.vcf.gz", false);
}

//with extended target region
$out_file3 = output_folder().$name."_out3.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file3 -model_type WES -target ".data_folder().$name."_in.bed -min_af 0.21 -target_extend 50 --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.vcf.gz");


########################## raw output ##########################

$out_file4 = output_folder().$name."_out4.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file4 -model_type WES -target ".data_folder().$name."_in.bed -raw_output --log ".output_folder().$name."_out4.log");
remove_lines_containing($out_file4, ["##contig=", "##fileDate=", "##commandline=", "##reference=", "11067308", "11065934"]); //the last two entries are genomic positions of variants the cause numeric problems 
check_file($out_file4, data_folder().$name."_out4.vcf.gz");

########################## gvcf output ##########################

$out_file6 = output_folder().$name."_out6.vcf.gz";
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in.bam -out $out_file6 -model_type WES -target ".data_folder().$name."_in.bed -gvcf ".output_folder().$name."_out6.gvcf.gz --log ".output_folder().$name."_out6.log");
check_file(output_folder().$name."_out6.gvcf.gz", data_folder().$name."_out6.gvcf.gz");

end_test();

?>


