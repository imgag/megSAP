<?php

require_once("framework.php");

$name = "an_vep";
start_test($name);

// function which sorts the Consequence entries of the VEP annotation
function sort_consequences($filename)
{
	$file_content = file($filename);
	for ($i=0; $i < count($file_content); $i++) 
	{
		if (substr($file_content[$i], 0, 14) === "chr12	11420421")
		{
			// get INFO column:
			$columns = explode("\t", $file_content[$i]);
			$info_column = $columns[7];

			//get VEP entries:
			$info_entries = explode(";", $info_column);
			for ($j=0; $j < count($info_entries); $j++) 
			{
				if (substr($info_entries[$j], 0, 11) == "CSQ_refseq=")
				{
					// get transcript
					$transcripts = explode(",", $info_entries[$j]);
					for ($k=0; $k < count($transcripts); $k++)
					{
						// get Consequence field:
						$csq = explode("|", $transcripts[$k]);
						$consequences = explode("&", $csq[1]);
					
						//sort
						sort($consequences);
						$csq[1] = implode("&", $consequences);

						$transcripts[$k] = implode("|", $csq);
					}
					

					// put line back together
					$info_entries[$j] = implode(",", $transcripts);

				}
			}
			$columns[7] = implode(";", $info_entries);
			$file_content[$i] = implode("\t", $columns);
		}
	}

	// write modified array to file
	file_put_contents($filename, $file_content);

}

// following tests will be run without group annotation which require NGSD connection (not available in nightly tests)

//standard
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in1.vcf -out $out_file1 --log ".output_folder().$name."_out1.log -no_groups");
remove_lines_containing($out_file1, array("##VEP=\"v"));

// TODO: Remove when VEP generates consistent refseq output for this variant
sort_consequences($out_file1);
//remove_lines_containing($out_file1, array("chr12	11420421	.	TCCTCCTTGTGGGGGTGGTCCTTCTGGCTTTCCTGGACGAGGTGGGGGACCTTGAGGTTTGTTG"));
check_file($out_file1, data_folder().$name."_out1.vcf", true);

//empty input VCF
$out_file2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -out $out_file2 --log ".output_folder().$name."_out2.log -no_groups");
remove_lines_containing($out_file2, array("##VEP=\"v"));
check_file($out_file2, data_folder().$name."_out2.vcf", true);

//somatic
$out_file3 = output_folder().$name."_out3.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -somatic -in ".data_folder().$name."_in1.vcf -out $out_file3 --log ".output_folder().$name."_out3.log -no_groups");
remove_lines_containing($out_file3, array("##VEP=\"v"));
// TODO: Remove when VEP generates consistent refseq output for this variant
sort_consequences($out_file3);
//remove_lines_containing($out_file3, array("chr12	11420421	.	TCCTCCTTGTGGGGGTGGTCCTTCTGGCTTTCCTGGACGAGGTGGGGGACCTTGAGGTTTGTTG"));
check_file($out_file3, data_folder().$name."_out3.vcf", true);

//NA12878_38 head
$out_file4 = output_folder().$name."_out4.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in3.vcf -out $out_file4 --log ".output_folder().$name."_out4.log -no_groups");
remove_lines_containing($out_file4, array("##VEP=\"v"));
check_file($out_file4, data_folder().$name."_out4.vcf", true);

//NA12878_38 head zipped
$out_file5 = output_folder().$name."_out5.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in4.vcf.gz -out $out_file5 --log ".output_folder().$name."_out5.log -no_groups");
remove_lines_containing($out_file5, array("##VEP=\"v"));
check_file($out_file5, data_folder().$name."_out5.vcf", true);

if (db_is_enabled("NGSD"))
{
	// these test will be run only if NGSD is available (not in nightly tests)

	//standard
	$out_file6 = output_folder().$name."_out6.vcf";
	check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in1.vcf -out $out_file6 --log ".output_folder().$name."_out6.log");
	remove_lines_containing($out_file6, array("##VEP=\"v"));
	// TODO: Remove when VEP generates consistent refseq output for this variant
	sort_consequences($out_file6);
	//remove_lines_containing($out_file6, array("chr12	11420421	.	TCCTCCTTGTGGGGGTGGTCCTTCTGGCTTTCCTGGACGAGGTGGGGGACCTTGAGGTTTGTTG"));
	check_file($out_file6, data_folder().$name."_out6.vcf", true);

	//NA12878_38 head
	$out_file7 = output_folder().$name."_out7.vcf";
	check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in3.vcf -out $out_file7 --log ".output_folder().$name."_out7.log");
	remove_lines_containing($out_file7, array("##VEP=\"v"));
	check_file($out_file7, data_folder().$name."_out7.vcf", true);

}

end_test();

?>
