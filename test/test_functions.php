<?php

require_once("framework.php");


//##################################################################################
start_test("container_platform");

check(container_platform(), "apptainer");
check(starts_with(container_platform(true), "apptainer 1."), true);

end_test();


//##################################################################################
start_test("is_in_temp_folder");

check(is_in_temp_folder("/tmp/local_ngs_data_GRCh38/"), true);
check(is_in_temp_folder("/tmp/non_existing_file"), true);
check(is_in_temp_folder("/home/"), false);

end_test();

//##################################################################################
start_test("chr2NC");

check(chr2NC("chr1"), "NC_000001.11");
check(chr2NC("chr11"), "NC_000011.10");
check(chr2NC("chrMT"), "NC_012920.1");

check(chr2NC("NC_000001.11", true), "chr1");
check(chr2NC("NC_000011.10", true), "chr11");
check(chr2NC("NC_012920.1", true), "chrMT");
end_test();

//##################################################################################
start_test("bound");

check(bound(-1.0,0.0,1.0), 0.0);
check(bound(0.0,0.0,1.0), 0.0);
check(bound(0.5,0.0,1.0), 0.5);
check(bound(1.0,0.0,1.0), 1.0);
check(bound(2.0,0.0,1.0), 1.0);

end_test();

//##################################################################################
start_test("allele_count");
check(allele_count(data_folder()."/allele_count_in.bam","chr1",31905889),array("C"=>0, "T"=>0, "A"=>263, "G"=>0, "*"=>110));
check(allele_count(data_folder()."/allele_count_in.bam","chr5",142421388),array("C"=>501, "T"=>336, "A"=>1, "G"=>0, "*"=>0));
check(allele_count(data_folder()."/allele_count_in.bam","chr13",110436093),array("C"=>0, "T"=>0, "A"=>319, "G"=>42, "*"=>1)); 
check(allele_count(data_folder()."/allele_count_in.bam","chr19",1221301),array("C"=>192, "T"=>0, "A"=>0, "G"=>4, "*"=>3));
end_test();
//##################################################################################
start_test("common_prefix");

check(common_prefix("abc", "xyz"), "");
check(common_prefix("", "xyz"), "");
check(common_prefix("abc", ""), "");
check(common_prefix("abc", "ab"), "ab");
check(common_prefix("abc", "abx"), "ab");
check(common_prefix("ab", "abc"), "ab");
check(common_prefix("abc", "abc"), "abc");

end_test();

//##################################################################################
start_test("common_suffix");

check(common_suffix("abc", "xyz"), "");
check(common_suffix("", "xyz"), "");
check(common_suffix("abc", ""), "");
check(common_suffix("abc", "bc"), "bc");
check(common_suffix("abc", "xbc"), "bc");
check(common_suffix("bc", "abc"), "bc");
check(common_suffix("abc", "abc"), "abc");

end_test();

//##################################################################################
start_test("repository_revision");

$rev = repository_revision();
check(starts_with($rev, "2025_10"), true);

$rev = repository_revision(true);
check(starts_with($rev, "megSAP 2025_10-"), true);

end_test();

//##################################################################################
start_test("time_readable");

check(time_readable(15.33333), "15.3333s");
check(time_readable(60), "1m 0s");
check(time_readable(115.33333), "1m 55s");
check(time_readable(3600), "1h 0m 0s");
check(time_readable(10015.33333), "2h 46m 55s");

end_test();

//##################################################################################
start_test("median");

check(median(array(4,1,8)), 4);
check(median(array(1,4,2,8)), 3);

end_test();


//##################################################################################
start_test("mean");

$data = array(1,1,1);
check(mean($data), 1);
$data = array(1,3,8);
check(mean($data), 4);

end_test();


//##################################################################################
start_test("stdev");

$data = array(1,1,1);
check(stdev($data), 0);
$data = array(3.2,3.5,2.9,3.3,3.4,2.5,2.7,2.8,3.1,2.6);
check(stdev($data), 0.331, 0.001);

end_test();

//##################################################################################
start_test("range_overlap");

check(range_overlap(-1, 5, 0, 3), true);
check(range_overlap(-1, 5, 4, 6), true);
check(range_overlap(-1, 5, -5, 0), true);
check(range_overlap(-1, 5, 6, 8), false);
check(range_overlap(-1, 5, -4, -2), false);

end_test();

//##################################################################################
start_test("random_string");

check(random_string(4, "A"), "AAAA");
check(random_string(2, "B"), "BB");

end_test();

//##################################################################################
start_test("get_timestamp");

$time1 = get_timestamp();
sleep(1);
$time2 = get_timestamp();

check(strlen($time1), 24);
check($time1!=$time2, true);

end_test();

//##################################################################################
start_test("starts_with");

check(starts_with("abc","a"), true);
check(starts_with("","a"), false);
check(starts_with("abc","b"), false);

end_test();

//##################################################################################
start_test("contains");

check(contains("abc","a"), true);
check(contains("","a"), false);
check(contains("abc","c"), true);
check(contains("abc","abc"), true);

end_test();

//##################################################################################
start_test("ends_with");

check(ends_with("abc","c"), true);
check(ends_with("","a"), false);
check(ends_with("abc","b"), false);

end_test();

//##################################################################################
start_test("nl_trim");

check(nl_trim("abc\n"), "abc");
check(nl_trim(" abc\t\n"), " abc\t");

end_test();

//##################################################################################
start_test("correlation");

$array1 = array("0.88", "1.04", "1.03", "1.10", "0.91", "0", "1.34", "1.15", "0.99", "1.14", "1.08");
$array2 = array("0.87", "1.05", "1.01", "1.14", "0.90", "0.02", "1.35", "1.10", "0.94", "1.17", "1.02");
$array3 = array("-0.87", "-1.05", "-1.01", "-1.14", "-0.90", "-0.02", "-1.35", "-1.10", "-0.94", "-1.17", "-1.02");
$array4 = array("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11");

check(correlation($array1, $array2), 1, 0.01);
check(correlation($array1, $array3), -1, 0.01);
check(correlation($array1, $array4), 0, 0.2);

end_test();

//##################################################################################
start_test("temp_file");

$file = temp_file(".suff", "pre");
$basename = basename($file);
$dir = dirname($file);

check(strlen($basename)>0, true);
check(strlen($dir)>0, true);
check(starts_with($basename, "pre"), true);
check(ends_with($basename, ".suff"), true);

end_test();

//##################################################################################
start_test("temp_folder");

$folder = temp_folder();

check(strlen($dir)>0, true);
check(is_dir($dir), true);
check(is_writable($dir), true);

end_test();

//##################################################################################
start_test("xml_is_wellformed");

check(xml_is_wellformed("<a><b></b><b><c></c></b></a>", $messages), true);
check(count($messages), 0);

check(xml_is_wellformed("<a><b></b><b><c></b></a>", $messages), false);
check(count($messages), 3);

end_test();

//##################################################################################
start_test("xml_matches_schema");

$schema =  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
			<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">
				<xs:element name=\"a\">
				  <xs:complexType>
					<xs:sequence>
					  <xs:element name=\"b\" maxOccurs=\"unbounded\"/>
					</xs:sequence>
				  </xs:complexType>
				</xs:element>
			</xs:schema>";

check(xml_matches_schema("<a><b></b><b><c></c></b></a>", $schema, $messages), true);
check(count($messages), 0);

check(xml_matches_schema("<a><b></b><b></b><c></c></a>", $schema, $messages), false);
check(count($messages), 1);

end_test();

//##################################################################################
start_test("xml_print_messages");

$lines = "<a><b></b><b><c></b></a>";
xml_is_wellformed($lines, $messages);
$output = xml_print_messages($messages, null, false);

check(count($output), 3);
check(starts_with($output[0], "Fatal"), true);
check(starts_with($output[1], "Fatal"), true);
check(starts_with($output[2], "Fatal"), true);

end_test();

//##################################################################################
start_test("range_intersect");

check(range_intersect(1, 2, 3, 4), false);
check(range_intersect(1, 2, -4, -2), false);
check(range_intersect(1, 5, 1, 5), array(1,5));
check(range_intersect(-1, 2, 1, 5), array(1,2));
check(range_intersect(1, 5, 4, 7), array(4,5));

end_test();


//##################################################################################
start_test("array_subset");

check(array_subset(array("A", "B", "C", "D"), 0, 1), array("A", "B"));
check(array_subset(array("A", "B", "C", "D"), 0, 3), array("A", "D"));
check(array_subset(array("A", "B", "C", "D"), 0, 0), array("A", "A"));
check(array_subset(array("A", "B", "C", "D"), 3, 2, 1, 0), array("D", "C", "B", "A"));
check(array_subset(array("A", "B", "C", "D"), array(0, 1)), array("A", "B"));
check(array_subset(array("A", "B", "C", "D"), array(0, 3)), array("A", "D"));
check(array_subset(array("A", "B", "C", "D"), array(0, 0)), array("A", "A"));
check(array_subset(array("A", "B", "C", "D"), array(3, 2, 1, 0)), array("D", "C", "B", "A"));

end_test();


//##################################################################################
start_test("mad");

$data = array(1,1,1);
check(stdev($data), 0);
$data = array(1, 1, 2, 2, 4, 6, 9);
check(mad($data), 1.4826, 0.01);
$data = array(-3, 15, 1, 1.4, 2.5, 2, 4.4, 6.5, 9.7, 4.2, 3.5, 5.5);
check(mad($data), 3.1876, 0.01);

end_test();


//##################################################################################
start_test("load_tsv");

$output = array();
$output[] = "#bla\n";
$output[] = "\n";
$output[] = "one\ttwo\n";
$output[] = " three\t four \tfive \n";
$output[] = "\n";
$output[] = "six\n";
$output[] = "#bluff\n";
file_put_contents("/tmp/php_test_load_tsv.txt", $output);

$data = load_tsv("/tmp/php_test_load_tsv.txt");

check(count($data), 3);
check($data[0], array("one", "two"));
check($data[1], array(" three", " four ", "five "));
check($data[2], array("six"));

unlink("/tmp/php_test_load_tsv.txt");

end_test();

//##################################################################################
start_test("sort_vcf_comments");

$comments = array(
	"#contig=<ID=chrUn_gl000249,length=38502>",
	"#SAMPLE=<ID=DX123458_01,Gender=n/a,DiseaseStatus=affected>",
	"#SAMPLE=<ID=DX123457_01,Gender=female,DiseaseStatus=control>",
	"#SAMPLE=<ID=DX123456_01,Gender=male,DiseaseStatus=control>",
	"#content=strelka somatic indel calls",
	"#INFO=<ID=QSI,Number=1,Type=Integer,Description=\"Quality score for any somatic variant, ie. for the ALT haplotype to be present at a significantly different frequency in the tumor and normal\">",
	"#SnpSiftCmd=\"SnpSift dbnsfp -f Interpro_domain /tmp/annotate_XwplOa_somatic.vcf\"",
	"#INFO=<ID=dbNSFP_Interpro_domain,Number=A,Type=String,Description=\"Field 'Interpro_domain' from dbNSFP\">",
	);
$comments_out = array(
	"#content=strelka somatic indel calls",
	"#SnpSiftCmd=\"SnpSift dbnsfp -f Interpro_domain /tmp/annotate_XwplOa_somatic.vcf\"",
	"#contig=<ID=chrUn_gl000249,length=38502>",
	"#INFO=<ID=dbNSFP_Interpro_domain,Number=A,Type=String,Description=\"Field 'Interpro_domain' from dbNSFP\">",
	"#INFO=<ID=QSI,Number=1,Type=Integer,Description=\"Quality score for any somatic variant, ie. for the ALT haplotype to be present at a significantly different frequency in the tumor and normal\">",
	"#SAMPLE=<ID=DX123458_01,Gender=n/a,DiseaseStatus=affected>",
	"#SAMPLE=<ID=DX123457_01,Gender=female,DiseaseStatus=control>",
	"#SAMPLE=<ID=DX123456_01,Gender=male,DiseaseStatus=control>",
);
check(sort_vcf_comments($comments), $comments_out);

end_test();


//##################################################################################
start_test("add_missing_contigs_to_vcf");

$input = data_folder()."/add_missing_contigs_to_vcf_in.vcf";
$out_grch38_ref = data_folder()."/add_missing_contigs_to_vcf_out_grch38.vcf";
$out_grch38 = output_folder()."/add_missing_contigs_to_vcf_out_grch38.vcf";

//copy files
copy($input, $out_grch38);

//annotate contigs
add_missing_contigs_to_vcf("GRCh38", $out_grch38);

// compare 
check_file($out_grch38, $out_grch38_ref, true);

end_test();

//##################################################################################
start_test("resolve_symlink");

//create symlink
$folder = realpath(data_folder());
check_exec("rm -rf {$folder}/symlink_file.txt {$folder}/real_file.txt"); //cleanup of previous tests
check_exec("touch {$folder}/real_file.txt");
check_exec("ln -s {$folder}/real_file.txt {$folder}/symlink_file.txt");

//perform checks
check(resolve_symlink("{$folder}/missing_file.txt"), "{$folder}/missing_file.txt");
check(resolve_symlink("{$folder}/real_file.txt"), "{$folder}/real_file.txt");
check(resolve_symlink("{$folder}/symlink_file.txt"), "{$folder}/real_file.txt");
check_exec("rm -rf {$folder}/symlink_file.txt {$folder}/real_file.txt"); //cleanup of previous tests
  
end_test();

?>
