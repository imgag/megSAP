<?php

require_once("framework.php");

$name = "an_vep";
start_test($name);
//TODO: Remove when VEP generates consistent refseq output for this variant
//function which sorts the Consequence entries of the VEP annotation
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

//standard
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in1.vcf -out $out_file1 --log ".output_folder().$name."_out1.log");
remove_lines_containing($out_file1, array("##VEP=\"v"));
remove_lines_containing($out_file1, array("##INFO=<ID=NGSD_GENE_INFO,"));
sort_consequences($out_file1);
check_file($out_file1, data_folder().$name."_out1.vcf", true);

//gzipped input VCF
$out_file2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in2.vcf.gz -out $out_file2 --log ".output_folder().$name."_out2.log");
remove_lines_containing($out_file2, array("##VEP=\"v"));
remove_lines_containing($out_file2, array("##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file2, data_folder().$name."_out2.vcf", true);

//empty input VCF
$out_file_empty = output_folder().$name."_out_empty.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -out $out_file_empty --log ".output_folder().$name."_out_empty.log");
remove_lines_containing($out_file_empty, array("##VEP=\"v"));
remove_lines_containing($out_file_empty, array("##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file_empty, data_folder().$name."_out_empty.vcf", true);

//variants to score with splicing variants
$out_file3 = output_folder().$name."_out3.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in3.vcf -out $out_file3 --log ".output_folder().$name."_out3.log");
remove_lines_containing($out_file3, array("##VEP=\"v"));
remove_lines_containing($out_file3, array("##INFO=<ID=NGSD_GENE_INFO,"));
sort_consequences($out_file3);
check_file($out_file3, data_folder().$name."_out3.vcf", true);

//DRAGEN
$out_file_dragen = output_folder().$name."_out_dragen.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in_dragen.vcf -out $out_file_dragen --log ".output_folder().$name."_out_dragen.log");
remove_lines_containing($out_file_dragen, array("##VEP=\"v"));
remove_lines_containing($out_file_dragen, array("##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file_dragen, data_folder().$name."_out_dragen.vcf", true);

//tests with NGSD disease group annotations - these tests will be run only if NGSD is available (i.e. not in nightly tests)
if (db_is_enabled("NGSD")) 
{
	//standard
	$out_file_db1 = output_folder().$name."_out_db1.vcf";
	check_exec("php ".src_folder()."/NGS/{$name}.php -test -ps_name NA12878_03 -in ".data_folder().$name."_in1.vcf -out $out_file_db1 --log ".output_folder().$name."_out_db1.log");
	remove_lines_containing($out_file_db1, array("##VEP=\"v"));
	remove_lines_containing($out_file_db1, array("##INFO=<ID=NGSD_GENE_INFO,"));
	sort_consequences($out_file_db1);
	check_file($out_file_db1, data_folder().$name."_out_db1.vcf", true);

	//somatic
	$out_file_som = output_folder().$name."_out_som.vcf";
	check_exec("php ".src_folder()."/NGS/{$name}.php -test -somatic -in ".data_folder().$name."_in1.vcf -out $out_file_som --log ".output_folder().$name."_out_som.log");
	remove_lines_containing($out_file_som, array("##VEP=\"v"));
	remove_lines_containing($out_file_som, array("##INFO=<ID=NGSD_GENE_INFO,"));
	sort_consequences($out_file_som);
	check_file($out_file_som, data_folder().$name."_out_som.vcf", true);
}

//Somatic VICC data
$out_file_som_vicc = output_folder().$name."_out_som_vicc.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -test -in ".data_folder().$name."_in5_som_vicc.vcf -out $out_file_som_vicc -somatic -no_splice --log ".output_folder().$name."_out_som_vicc.log");
remove_lines_containing($out_file_som_vicc, array("##VEP=\"v", "##INFO=<ID=NGSD_GROUP"));
remove_lines_containing($out_file_som_vicc, array("##INFO=<ID=NGSD_GENE_INFO,"));
check_file($out_file_som_vicc, data_folder().$name."_out_som_vicc.vcf", true);

end_test();

?>
