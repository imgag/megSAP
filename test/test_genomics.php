<?php

include("framework.php");


//##################################################################################
start_test("bed_is_sorted");

check(bed_is_sorted(data_folder()."/an_vep_NGSD_gene_info.bed"), true);
check(bed_is_sorted(data_folder()."/bed_unsorted.bed"), false);

end_test();

//##################################################################################
start_test("get_qcml_value");

$qcml = data_folder()."/test.qcML";

check(get_qcml_value($qcml, "QC:2000025"), "180.11");
check(get_qcml_value($qcml, "QC:2000027"), "93.46");

check(get_qcml_value($qcml, "QC:2000026", "n/a"), "n/a");

end_test();

//##################################################################################
start_test("genome_fasta");

$genome_local = genome_fasta("GRCh38");
check(file_exists($genome_local), true); //data from megSAP folder
$genome_repo = genome_fasta("GRCh38", false);
check(file_exists($genome_repo), true); //local copy of data
check($genome_local!=$genome_repo, true);

end_test();

//##################################################################################
start_test("vcfgeno2human");

check(vcfgeno2human("0|0"), "wt");
check(vcfgeno2human("./1"), "het");
check(vcfgeno2human("1|1"), "hom");
check(vcfgeno2human("1|1", true), "HOM");

end_test();

//##################################################################################
start_test("is_valid_ref_sample_for_cnv_analysis");

check(is_valid_ref_sample_for_cnv_analysis("NA12878"), false);
if (db_is_enabled("NGSD") && get_path("location", false)=="IMGAG")
{
	check(is_valid_ref_sample_for_cnv_analysis("GS160265_06"), false); //tumor
	check(is_valid_ref_sample_for_cnv_analysis("GS160832_01"), false); //ffpe
	check(is_valid_ref_sample_for_cnv_analysis("I16D005a01_01"), false); //bad run
	check(is_valid_ref_sample_for_cnv_analysis("GS130071_01"), false); //bad processed sample
	check(is_valid_ref_sample_for_cnv_analysis("GS123456_01"), false); //not in NGSD
	check(is_valid_ref_sample_for_cnv_analysis("GS150155_05"), false); //test project
	check(is_valid_ref_sample_for_cnv_analysis("GS150155_05", false, true), true); //test project
	check(is_valid_ref_sample_for_cnv_analysis("GS160408_01"), true);
}
end_test();

//##################################################################################
start_test("correct_indel");

//no change
check(correct_indel(5,"A","G"), array(5,5,"A","G"));
check(correct_indel(5,"A","CG"), array(5,5,"A","CG"));

//start del + ins
check(correct_indel(5,"AcG","CG"), array(5,5,"A","-"));
check(correct_indel(5,"CG","AcG"), array(4,4,"-","A"));
check(correct_indel(5,"ACGCG","CGCG"), array(5,5,"A","-"));
check(correct_indel(5,"CGCG","ACGCG"), array(4,4,"-","A"));

//middle del + ins
check(correct_indel(5,"ACG","AG"), array(6,6,"C","-"));
check(correct_indel(5,"AG","ACG"), array(5,5,"-","C"));
check(correct_indel(5,"AACGG","AAGG"), array(7,7,"C","-"));
check(correct_indel(5,"AAGG","AACGG"), array(6,6,"-","C"));
check(correct_indel(5,"AAACG","AAAG"), array(8,8,"C","-"));
check(correct_indel(5,"AAAG","AAACG"), array(7,7,"-","C"));
check(correct_indel(5,"ACGGG","AGGG"), array(6,6,"C","-"));
check(correct_indel(5,"AGGG","ACGGG"), array(5,5,"-","C"));

//end del + ins
check(correct_indel(5,"ACG","AC"), array(7,7,"G","-"));
check(correct_indel(5,"AC","ACG"), array(6,6,"-","G"));
check(correct_indel(5,"ACACG","ACAC"), array(9,9,"G","-"));
check(correct_indel(5,"ACAC","ACACG"), array(8,8,"-","G"));

//complex indels
check(correct_indel(5,"ACGT","TGCA"), array(5,8,"ACGT","TGCA"));
check(correct_indel(2832467,"GGGGGCG","C"), array(2832467,2832473,"GGGGGCG","C"));

end_test();

//##################################################################################
start_test("gender");

$geno_f = array("AB", "BB", "AA", "BB", "AB", "AA", "AB", "BB", "AB", "BB", "AA", "BB", "AB", "AA", "AB", "BB");
$geno_m = array("AA", "BB", "AA", "BB", "AB", "AA", "BB", "BB", "AA", "BB", "AA", "BB", "AA", "AA", "AB", "BB");

check(gender($geno_f, "AB", 0.3, 0.3), array('f', 0.375));
check(gender($geno_m, "AB", 0.3, 0.3), array('m', 0.125));
check(gender($geno_f, "AB", 0.3, 1.0), array(false, 0.375));
check(gender($geno_m, "AB", 0.0, 0.3), array(false, 0.125));
check(gender($geno_f, "AB", 0.7, 0.7), array('m', 0.375));
check(gender($geno_m, "AB", 0.0, 0.0), array('f', 0.125));
end_test();

//##################################################################################
start_test("rev_comp");

check(rev_comp("ACGT"), "ACGT");
check(rev_comp("AA"), "TT");
check(rev_comp("ac"), "gt");

end_test();

//##################################################################################
start_test("from_IUPAC");

check(from_IUPAC("R"), array("A", "G"));

end_test();

//##################################################################################
start_test("chr_trim");

check(chr_trim("X"), "X");
check(chr_trim("x"), "X");
check(chr_trim("Y"), "Y");
check(chr_trim("y"), "Y");
check(chr_trim("chrX"), "X");
check(chr_trim("chrM"), "M");
check(chr_trim("chr1"), "1");
check(chr_trim("CHR1"), "1");
check(chr_trim("1"), "1");
check(chr_trim("15"), "15");
check(chr_trim("22"), "22");
check(chr_trim("chr22"), "22");
check(chr_trim("CHR22"), "22");

end_test();

//##################################################################################
start_test("chr_check");

check(chr_check("x"), true);
check(chr_check("Y"), true);
check(chr_check("chrM"), true);
check(chr_check("1"), true);
check(chr_check("chr22"), true);

check(chr_check("chr0",22,false), false);
check(chr_check("0",22,false), false);
check(chr_check("chr23",22,false), false);
check(chr_check("23",22,false), false);
check(chr_check("chrO",22,false), false);

end_test();

//##################################################################################
start_test("get_ref_seq");

check(get_ref_seq("GRCh38", "chr14", 57349540, 57349542), "AAT");
check(get_ref_seq("GRCh38", "chr14", 57349543, 57349545), "GCC");
check(get_ref_seq("GRCh38", "chr12", 12310000, 12310200), "CACCATGCCCGGCTGGTGGAGTTTTTTTAACCTAAAAATTTTACTTCTGTATTAACTTTATAGCACAAATATTTTCTTCATAATACGGATAATCCAGAGGGAAGTTCAAAAGTGAAACAGATGGAGGAAGGGCAGGGGAAAAGATGGTCCCTTGAATTGGGGAGAGTAGAAGGTCTCAGGGTCCCCATGTCATGGCTTCCC");

// test caching
check(get_ref_seq("GRCh38", "chr14", 57349540, 57349542, 100), "AAT");
check(get_ref_seq("GRCh38", "chr14", 57349543, 57349545, 100), "GCC");
check(get_ref_seq("GRCh38", "chr12", 12310000, 12310200, 100), "CACCATGCCCGGCTGGTGGAGTTTTTTTAACCTAAAAATTTTACTTCTGTATTAACTTTATAGCACAAATATTTTCTTCATAATACGGATAATCCAGAGGGAAGTTCAAAAGTGAAACAGATGGAGGAAGGGCAGGGGAAAAGATGGTCCCTTGAATTGGGGAGAGTAGAAGGTCTCAGGGTCCCCATGTCATGGCTTCCC");

// test behavior at chr ends
check(get_ref_seq("GRCh38", "chr1", 1, 10, 100), "NNNNNNNNNN");
check(get_ref_seq("GRCh38", "chr1", 21, 30, 100), "NNNNNNNNNN");
check(get_ref_seq("GRCh38", "chr1", 248956422, 248956422, 1000), "N");
check(get_ref_seq("GRCh38", "chr1", 248956000, 248956000, 1000), "N");
check(get_ref_seq("GRCh38", "chrX", 156030895, 156030895, 100000), "T");
check(get_ref_seq("GRCh38", "chrX", 156030895, 156030895, 100000), "T");
check(get_ref_seq("GRCh38", "chrX", 156030895, 156030895, 10), "T");
check(get_ref_seq("GRCh38", "chrX", 156030895, 156030895, 0), "T");

// test on non-local genome file
check(get_ref_seq("GRCh38", "chr14", 57349540, 57349542, 0, false), "AAT");

end_test();

//##################################################################################
start_test("load_system");

$filename = data_folder()."/genomics_load_system.ini";
$sys = load_system($filename);
check($sys['name_manufacturer'], "SureSelectXT Human All Exon v5");
check($sys['shotgun'], true);
check($sys['umi_type'], "n/a");

if (db_is_enabled("NGSD") && get_path("location", false)=="IMGAG")
{
	$filename = "";
	$sys = load_system($filename, "GS130043_01");
	check($filename!="", true);
	check($sys['name_short'], "hpSCAv2");
}

end_test();

//##################################################################################
start_test("store_system");

if (db_is_enabled("NGSD") && get_path("location", false)=="IMGAG")
{
	$db_conn = DB::getInstance("NGSD");
	$filename = temp_file(".ini");
	$sys = store_system($db_conn, "hpSCAv2", $filename);
	check(file_exists($filename), true);
	check(contains(file_get_contents($filename), 'name_short = "hpSCAv2"'), true);
}

end_test();

//##################################################################################
start_test("get_processed_sample_id");

if (db_is_enabled("NGSD") && get_path("location", false)=="IMGAG")
{
	$db_conn = DB::getInstance("NGSD");
	check(get_processed_sample_id($db_conn, "GS130043_01"), 1498);
	check(get_processed_sample_id($db_conn, "GS130043_01_tumor"), 1498);
	check(get_processed_sample_id($db_conn, "GS130043", false), -1);
	check(get_processed_sample_id($db_conn, "GS123456_01", false), -1);
}

end_test();

//##################################################################################
start_test("gsvar_to_vcf");
check(gsvar_to_vcf("GRCh38", "chr12",25368401,"T","A"),array("chr12",25368401,"T","A"));
check(gsvar_to_vcf("GRCh38", "chr12",25368401,"TC","AG"),array("chr12",25368401,"TC","AG"));
check(gsvar_to_vcf("GRCh38", "chr12",25368401,"TC","-"),array("chr12",25368400,"TTC","T"));
check(gsvar_to_vcf("GRCh38", "chr12",25368401,"-","AG"),array("chr12",25368401,"G","GAG"));
check(gsvar_to_vcf("GRCh38", "chr18",50207290,"A","GTGTGTG"),array("chr18",50207289,"TA","TGTGTGTG"));	
end_test();

//##################################################################################
start_test("vcf_strelka_snv");
check(vcf_strelka_snv("AU:CU:DP:FDP:GU:SDP:SUBDP:TU","0,0:0,0:236:3:231,314:0:0:2,6","A"), array(236,0), 1e-3);
check(vcf_strelka_snv("AU:CU:DP:FDP:GU:SDP:SUBDP:TU","0,0:0,0:236:3:231,314:0:0:2,6","C"), array(236,0), 1e-3);
check(vcf_strelka_snv("AU:CU:DP:FDP:GU:SDP:SUBDP:TU","0,0:0,0:236:3:231,314:0:0:2,6","T"), array(236,0.0086), 1e-3);
check(vcf_strelka_snv("AU:CU:DP:FDP:GU:SDP:SUBDP:TU","0,0:0,0:236:3:231,314:0:0:2,6","G"), array(236,0.9914), 1e-3);
check(vcf_strelka_snv("AU:CU:DP:FDP:GU:SDP:SUBDP:TU","0,0:0,0:236:3:0,0:0:0:0,0","G"), array(236,null), 1e-3);
end_test();

//##################################################################################
start_test("vcf_strelka_snv_postcall");
check(vcf_strelka_snv_postcall("DP:AU:GU:CU:TU","100:50,50:25,25:25,25:0,0", "A", "G", 0.05), [ [ "C", 0.25 ] ], null);
check(vcf_strelka_snv_postcall("DP:AU:GU:CU:TU","100:25,25:25,25:25,25:25,25", "A", "G", 0.05), [ [ "T", 0.25 ], [ "C", 0.25 ] ], null);
end_test();

//##################################################################################
start_test("vcf_strelka_indel");
check(vcf_strelka_indel("DP:DP2:DP50:FDP50:SUBDP50:TAR:TIR:TOR","117:117:115.65:1.38:0.00:17,19:29,31:74,77"), array(117,0.6304), 1e-3);
check(vcf_strelka_indel("DP:DP2:DP50:FDP50:SUBDP50:TAR:TIR:TOR","275:275:255.49:9.44:0.00:372,401:0,0:16,3"), array(275,0.0000), 1e-3);
check(vcf_strelka_indel("DP:DP2:DP50:FDP50:SUBDP50:TAR:TIR:TOR","275:275:255.49:9.44:0.00:0,0:0,0:16,3"), array(275,null), 1e-3);
end_test();

//##################################################################################
start_test("vcf_freebayes");
check(vcf_freebayes("GT:GL:DP:RO:QR:AO:QA","0/1:-23.0878,0,-65.9409:29:21:752:8:276"), array(29,0.2759));
check(vcf_freebayes("GT:GL:DP:RO:QR:AO:QA","0/1:-23.0878,0,-65.9409:0:21:752:8:276"), array(0,null));
check(vcf_freebayes("GT:DP:AD:RO:QR:AO:QA:GL","0/0:9:8,1:8:305:1:16:0,-1.10929,-26.2088"), array(9,0.1111));
check(vcf_freebayes("GT:DP:AD:RO:QR:AO:QA:GL","0/1:14:10,4:10:327:4:132:-7.98514,0,-25.5222"), array(14,0.2857));
end_test();

//##################################################################################
start_test("vcf_varscan2");
check(vcf_varscan2("GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR","0/1:255:1560:1560:1070:489:31.37%:2.7462E-165:36:35:707:363:355:134"), array(1560, 0.3137));
check(vcf_varscan2("GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR", "0/1:52:200:200:183:17:8.5%:5.3509E-6:36:36:138:45:14:3"), array(200,0.085));
end_test();


//##################################################################################
start_test("vcf_iontorrent");
check(vcf_iontorrent("GT:GQ:DP:FDP:RO:FRO:AO:FAO:AF:SAR:SAF:SRF:SRR:FSAR:FSAF:FSRF:FSRR","6/6:125:5162:3914:72:0:0,5,0,6,2,5072:3,2,0,0,0,3909:0.000766479,0.000510986,0,0,0,0.998723:0,3,0,3,0,2683:0,2,0,3,2,2389:35:37:0,0,0,0,0,2045:3,2,0,0,0,1864:0:0",5), array(5162,0.9987));
end_test();

//##################################################################################
start_test("vcf_column_index");
check(vcf_column_index("GS140549", array("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GS140549", "GS140127", "GS140550")), 9);
check(vcf_column_index("GS140127", array("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GS140549", "GS140127", "GS140550")), 10);
check(vcf_column_index("GS140550", array("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GS140549", "GS140127", "GS140550")), 11);
check(vcf_column_index("GS140549_01", array("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GS140549", "GS140127", "GS140550")), 9);
check(vcf_column_index("GS140127_01", array("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GS140549", "GS140127", "GS140550")), 10);
check(vcf_column_index("GS140550_01", array("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GS140549", "GS140127", "GS140550")), 11);
check(vcf_column_index("KLO0051", array("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "KLO0052_01", "KLO0051_01")), 10);
check(vcf_column_index("KLO0051_01", array("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "KLO0052", "KLO0051")), 10);
end_test();

start_test("aa1_to_aa3");
check( aa1_to_aa3("*") , "Ter" );
check( aa1_to_aa3("F") , "Phe" );
check( aa3_to_aa1("Ter") , "*" );
check( aa3_to_aa1("Asp") , "D" );
end_test();

//##################################################################################
start_test("get_genome_build");
check(check_genome_build(data_folder()."/get_genome_build_bwaGRCh37.bam", "GRCh37"), 1);
check(check_genome_build(data_folder()."/get_genome_build_bwaGRCh38.bam", "GRCh38"), 1);
check(check_genome_build(data_folder()."/get_genome_build_bwamem2GRCh37.bam", "GRCh37"), 1);
check(check_genome_build(data_folder()."/get_genome_build_bwamem2GRCh38.bam", "GRCh38"), 1);
check(check_genome_build(data_folder()."/get_genome_build_dragenGRCh37.bam", "GRCh37"), 1);
check(check_genome_build(data_folder()."/get_genome_build_dragenGRCh38.bam", "GRCh38"), 1);
check(check_genome_build(data_folder()."/get_genome_build_NovaSeqXGRCh38.bam", "GRCh38"), 1);
check(check_genome_build(data_folder()."/check_genome_build_minimap2GRCh38.bam", "GRCh38"), 1);
check(check_genome_build(data_folder()."/check_genome_build_minimap2GRCh38_2.bam", "GRCh38"), 1);
end_test();

//##################################################################################
start_test("get_processed_sample_info");
if (db_is_enabled("NGSD_TEST"))
{
	$db_conn = DB::getInstance("NGSD_TEST");
	$sql_file = data_folder()."/merge_samples.sql";
	execApptainer("ngs-bits", "NGSDInit", "-test -add $sql_file", [$sql_file]);
	
	$sample_info = get_processed_sample_info($db_conn, "DNA220002_01");
	check($sample_info["sys_target"], "");
	check($sample_info["ps_lanes"], array(1));
	check(ends_with($sample_info['project_folder'], "/merge_samples/"), true);
	check($sample_info['ps_name'], "DNA220002_01");
	check(ends_with($sample_info['ps_folder'], "/merge_samples/Sample_DNA220002_01/"), true);
	check(ends_with($sample_info['ps_bam'], "/merge_samples/Sample_DNA220002_01/DNA220002_01.bam"), true);

	$sample_info = get_processed_sample_info($db_conn, "DNA220002_01", true, true);
	check($sample_info["sys_target"], "");
	check($sample_info["ps_lanes"], array(1));
	check(isset($sample_info['project_folder']), false);
	check(isset($sample_info['ps_folder']), false);
	check(isset($sample_info['ps_bam']), false);
}
else
{
	print("No NGSD_TEST connection, skipping test!");
}
end_test();

//##################################################################################
start_test("basename2");

check(basename2("/some/path/filename1.bam")  , "filename1");
check(basename2("/some/path/filename2.BAM")  , "filename2");
check(basename2("/some/path/filename3.CRAM") , "filename3");
check(basename2("/some/path/filename4.cram") , "filename4");
check(basename2("/some/path/filename5")      , "filename5");

end_test();

//##################################################################################
start_test("get_basecall_model");

check(get_basecall_model(data_folder()."/get_basecall_model_dorado352.bam"), "dna_r10.4.1_e8.2_400bps_hac@v3.5.2");
check(get_basecall_model(data_folder()."/get_basecall_model_dorado420.bam"), "dna_r10.4.1_e8.2_400bps_hac@v4.2.0");
check(get_basecall_model(data_folder()."/get_basecall_model_no_model.bam"), "");

end_test();

//##################################################################################
start_test("chr_list");

check(chr_list(), array("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
						"chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
						"chr22", "chrX", "chrY"));

end_test();

//##################################################################################
start_test("get_read_counts");

check(get_read_count(data_folder()."/get_read_count_in1.bam"), 1861);
check(get_read_count(data_folder()."/get_read_count_in1.bam", 1, array("-F", "2304")), 1607);
check(get_read_count(data_folder()."/get_read_count_in1.bam", 1, array(), "GRCh38", "chr1:37999900-38000000"), 70);

check(get_read_count(data_folder()."/get_read_count_in2.cram"), 1861);
check(get_read_count(data_folder()."/get_read_count_in2.cram", 1, array("-F", "2304")), 1607);
check(get_read_count(data_folder()."/get_read_count_in2.cram", 1, array(), "GRCh38", "chr1:37988000-37999900"), 281);

end_test();

//##################################################################################
start_test("compare_bam_read_count");

check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in1.bam", data_folder()."/compare_bam_read_count_in1.cram"), true);
check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in1.cram", data_folder()."/compare_bam_read_count_in1.bam"), true);
check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in1.bam", data_folder()."/compare_bam_read_count_in2.bam"), true);
check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in1.bam", data_folder()."/compare_bam_read_count_in2.bam", 2, true, false, 0.0, array("-F", "2304")), true);
check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in1.bam", data_folder()."/compare_bam_read_count_in2.bam", 2, false, true), false);
check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in1.bam", data_folder()."/compare_bam_read_count_in2.bam", 2, false, true, 0.0, array("-F", "2304")), false);
check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in1.cram", data_folder()."/compare_bam_read_count_in2.cram", 2, false, true, 0.021), true);
check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in2.bam", data_folder()."/compare_bam_read_count_in3.bam", 2, false, false, 0.005), true);
check(compare_bam_read_count(data_folder()."/compare_bam_read_count_in2.bam", data_folder()."/compare_bam_read_count_in3.bam", 2, false, true, 0.041), true);

end_test();

//##################################################################################
start_test("get_ora_read_count");

check(get_ora_read_count([data_folder()."/get_ora_read_count_in1_L001_R1_001.fastq.ora"]), 1250);
check(get_ora_read_count([data_folder()."/get_ora_read_count_in1_L001_R2_001.fastq.ora"]), 1250);
check(get_ora_read_count([data_folder()."/get_ora_read_count_in1_L001_R1_001.fastq.ora", data_folder()."/get_ora_read_count_in1_L001_R2_001.fastq.ora"]), 2500);

end_test();

//##################################################################################
start_test("get_fastq_read_count");

check(get_fastq_read_count([data_folder()."/get_fastq_read_count_in1_L001_R1_001.fastq.gz"]), 1250);
check(get_fastq_read_count([data_folder()."/get_fastq_read_count_in1_L001_R2_001.fastq.gz"]), 1250);
check(get_fastq_read_count([data_folder()."/get_fastq_read_count_in1_L001_R1_001.fastq.gz", data_folder()."/get_fastq_read_count_in1_L001_R2_001.fastq.gz"]), 2500);

end_test();

//##################################################################################
start_test("contains_mito");

check(contains_mito(data_folder()."/an_vep_in1.vcf"), true);
check(contains_mito(data_folder()."/an_vep_in2.vcf"), false);

check(contains_mito(data_folder()."/vc_clair_out3.vcf.gz"), true);
check(contains_mito(data_folder()."/vc_clair_out2.vcf.gz"), false);

check(contains_mito(data_folder()."/an_vep_NGSD_gene_info.bed"), true);
check(contains_mito(data_folder()."/vc_clair_in_roi.bed"), false);

end_test();

?>
