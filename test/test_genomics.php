<?php

include("framework.php");

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
if (db_is_enabled("NGSD"))
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

check(get_ref_seq("GRCh37", "chr14", 57349540, 57349542), "ACT");
check(get_ref_seq("GRCh37", "chr14", 57349543, 57349545), "ATC");
check(get_ref_seq("GRCh37", "chr12", 12310000, 12310200), "AAGCCTATGAGAGAAAGCTGCTGGCTCTTGAACTATACCTTCTCTTTAGGTAACCTCATTCATTTTAAATACATCCTGGTAATCCCAAAATTTGTATCTTCAACCTCATGTCTCTTCTCCGATGATAGTCCATCATCTATGACACATAATCTAGAGGTCATCCTTGCTTCCTCCCTTTCTTTCCCACACAGAACTCATTAA");

// test caching
check(get_ref_seq("GRCh37", "chr14", 57349540, 57349542, 100), "ACT");
check(get_ref_seq("GRCh37", "chr14", 57349543, 57349545, 100), "ATC");
check(get_ref_seq("GRCh37", "chr12", 12310000, 12310200, 100), "AAGCCTATGAGAGAAAGCTGCTGGCTCTTGAACTATACCTTCTCTTTAGGTAACCTCATTCATTTTAAATACATCCTGGTAATCCCAAAATTTGTATCTTCAACCTCATGTCTCTTCTCCGATGATAGTCCATCATCTATGACACATAATCTAGAGGTCATCCTTGCTTCCTCCCTTTCTTTCCCACACAGAACTCATTAA");

// test behavior at chr ends
check(get_ref_seq("GRCh37", "chr1", 1, 10, 100), "NNNNNNNNNN");
check(get_ref_seq("GRCh37", "chr1", 21, 30, 100), "NNNNNNNNNN");
check(get_ref_seq("GRCh37", "chr1", 249250621, 249250621, 1000), "N");
check(get_ref_seq("GRCh37", "chr1", 249250021, 249250021, 1000), "N");
check(get_ref_seq("GRCh37", "chrX", 155260465, 155260465, 100000), "T");
check(get_ref_seq("GRCh37", "chrX", 155260465, 155260465, 100000), "T");
check(get_ref_seq("GRCh37", "chrX", 155260465, 155260465, 10), "T");
check(get_ref_seq("GRCh37", "chrX", 155260465, 155260465, 0), "T");

end_test();

//##################################################################################
start_test("load_system");

$filename = data_folder()."/genomics_load_system.ini";
$sys = load_system($filename);
check($sys['name_manufacturer'], "SureSelectXT Human All Exon v5");
check($sys['shotgun'], true);
check($sys['umi_type'], "n/a");

if (db_is_enabled("NGSD"))
{
	$filename = "";
	$sys = load_system($filename, "GS130043_01");
	check($filename!="", true);
	check($sys['name_short'], "hpSCAv2");
}

end_test();

//##################################################################################
start_test("store_system");

if (db_is_enabled("NGSD"))
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

if (db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	check(get_processed_sample_id($db_conn, "GS130043_01"), 1498);
	check(get_processed_sample_id($db_conn, "GS130043", false), -1);
	check(get_processed_sample_id($db_conn, "GS123456_01", false), -1);
}

end_test();

//##################################################################################
start_test("convert_coding2genomic");
#chr3	195508046	195508046	C	G	MUC4:NM_018406.6:missense:MODERATE:exon2/25:c.10405G>C:p.Asp3469His,MUC4:NM_004532.5:intron:MODIFIER:exon1/23:c.83-2720G>C:,MUC4:NM_138297.4:intron:MODIFIER:exon1/22:c.83-6870G>C:
#chr11	1092699	1092699	G	A	het	QUAL=7609;DP=1041;AF=0.27;MQM=56	MUC2	synonymous	MUC2:NM_002457.3:synonymous:LOW:exon30/50:c.4518G>A:p.Thr1506Thr
check(convert_coding2genomic("NM_018406.6",10405,10405),array("chr3",195508046,195508046,"-"));
check(convert_coding2genomic("NM_002457.3",4518,4518),array("chr11",1092699,1092699,"+"));
end_test();


//##################################################################################
start_test("convert_hgvs2genomic");

//SNVs
#--strand,coding and splicing	chr3	195508046	195508046	C	G	MUC4:NM_018406.6:missense:MODERATE:exon2/25:c.10405G>C:p.Asp3469His,MUC4:NM_004532.5:intron:MODIFIER:exon1/23:c.83-2720G>C:,MUC4:NM_138297.4:intron:MODIFIER:exon1/22:c.83-6870G>C:
check(convert_hgvs2genomic("GRCh37", "NM_018406.6", "c.10405G>C"),array("chr3",195508046,195508046,"C","G"));
check(convert_hgvs2genomic("GRCh37", "NM_004532.5", "c.83-2720G>C"),array("chr3",195508046,195508046,"C","G"));
#+-strand, coding	chr8	30924557	30924557	C	T	het	QUAL=7686;DP=561;AF=0.44;MQM=60	WRN	synonymous	WRN:NM_000553.4:synonymous:LOW:exon6/35:c.513C>T:p.Cys171Cys
check(convert_hgvs2genomic("GRCh37", "NM_000553.4", "c.513C>T"),array("chr8",30924557,30924557,"C","T"));
#--strand, coding	chr7	6013049	6013049	C	G	hom	QUAL=2763;DP=92;AF=1.00;MQM=35	PMS2	missense	PMS2:NM_000535.5:missense:MODERATE:exon15/15:c.2570G>C:p.Gly857Ala
check(convert_hgvs2genomic("GRCh37", "NM_000535.5", "c.2570G>C"),array("chr7",6013049,6013049,"C","G"));
#+-strand, splicing	chr9	17135434	17135434	G	C	het	QUAL=4662;DP=350;AF=0.44;MQM=60	CNTLN	splice_region&intron,intron	CNTLN:NM_017738.3:splice_region&intron:LOW:exon1/25:c.360+11G>C:,CNTLN:NM_001286984.1:splice_region&intron:LOW:exon1/2:c.360+11G>C:,CNTLN:NM_001114395.2:splice_region&intron:LOW:exon1/6:c.360+11G>C:,CNTLN:NM_001286985.1:intron:MODIFIER:exon1/2:c.171+200G>C:
check(convert_hgvs2genomic("GRCh37", "NM_017738.3", "c.360+11G>C"),array("chr9",17135434,17135434,"G","C"));
#+-strand, splicing	chr2	47630550	47630550	C	G	het	QUAL=2300;DP=168;AF=0.46;MQM=60	MSH2	splice_region&intron	MSH2:NM_000251.2:splice_region&intron:LOW:exon1/15:c.211+9C>G:,MSH2:NM_001258281.1:splice_region&intron:LOW:exon2/16:c.13+9C>G:
check(convert_hgvs2genomic("GRCh37", "NM_000251.2", "c.211+9C>G"),array("chr2",47630550,47630550,"C","G"));
#--strand;splicing:	chr18	21111585	21111585	C	T	het	QUAL=0	C18orf8,NPC1	splice_region&intron,3'UTR	C18orf8:NM_013326.4:splice_region&intron:LOW:exon19/19:c.1895-4C>T:,C18orf8:NM_001276342.1:splice_region&intron:LOW:exon17/17:c.1627-4C>T:,C18orf8:NR_075075.1:splice_region&intron:LOW:exon16/16:n.1815-4C>T:,C18orf8:NR_075076.1:splice_region&intron:LOW:exon18/18:n.1835-4C>T:,NPC1:NM_000271.4:3'UTR:MODIFIER:exon25/25:c.*581G>A:
check(convert_hgvs2genomic("GRCh37", "NM_013326.4", "c.1895-4C>T"),array("chr18",21111585,21111585,"C","T"));
#--strand,splicing:	chr18	21112158	21112158	C	A	NPC1	c.3838+8G>T
check(convert_hgvs2genomic("GRCh37", "NM_000271.4", "c.3837+8G>T"),array("chr18",21112158,21112158,"C","A"));


//Ins
#+-strand, splicing	chr11	108121410	108121410	-	T	het	QUAL=1626;DP=203;AF=0.30;MQM=57	ATM	splice_region&intron	ATM:NM_000051.3:splice_region&intron:LOW:exon9/62:c.1236-18_1236-17insT:
check(convert_hgvs2genomic("GRCh37", "NM_000051.3", "c.1236-18_1236-17insT"),array("chr11",108121410,108121410,"-","T"));
#c.3755-5_3755-4insTC, --strand,splicing:	chr18	21112158	21112158	C	A	NPC1	c.3838+8G>T
check(convert_hgvs2genomic("GRCh37", "NM_000271.4", "c.3755-5_3755-4insTC"),array("chr18",21112252,21112252,"-","GA"));
#--strand, splicing: chr3	142241692	142241692	-	A	het	QUAL=390;DP=391;AF=0.12;MQM=60	ATR	splice_region&intron	ATR:NM_001184.3:splice_region&intron:LOW:exon22/46:c.4153-10_4153-9insT:
check(convert_hgvs2genomic("GRCh37", "NM_001184.3", "c.4153-10_4153-9insT"),array("chr3",142241692,142241692,"-","A"));
#+-strand chr5	72743299	72743299	-	GC	hom	QUAL=206;DP=10;AF=0.80;MQM=60	FOXD1	frameshift	FOXD1:NM_004472.2:frameshift:HIGH:exon1/1:c.888_889insGC:p.Arg297fs
check(convert_hgvs2genomic("GRCh37", "NM_004472.2", "c.888_889insGC"),array("chr5",72743299,72743299,"-","GC"));
#chr3	75790810	75790810	-	T	het	QUAL=377;DP=20;AF=0.70;MQM=21	ZNF717,MIR4273	frameshift,5'UTR,downstream_gene	ZNF717:NM_001128223.1:frameshift:HIGH:exon3/5:c.134_135insA:p.Leu46fs,ZNF717:NM_001290210.1:frameshift:HIGH:exon3/6:c.134_135insA:p.Leu46fs,ZNF717:NM_001290208.1:frameshift:HIGH:exon3/5:c.134_135insA:p.Leu46fs,ZNF717:NM_001290209.1:5'UTR:MODIFIER:exon3/5:c.-17_-16insA:,MIR4273:NR_036235.1:downstream_gene:MODIFIER::n.*84_*84insT:
check(convert_hgvs2genomic("GRCh37", "NM_001128223.1", "c.134_135insA"),array("chr3",75790810,75790810,"-","T"));


//Del
#--strand, splicing	chr9	32986031	32986031	A	-	het	QUAL=2196;DP=105;AF=0.49;MQM=58	APTX	splice_region&intron	APTX:NM_001195248.1:splice_region&intron:LOW:exon4/7:c.526-3delT:,APTX:NM_001195249.1:splice_region&intron:LOW:exon4/7:c.484-3delT:,APTX:NM_001195254.1:splice_region&intron:LOW:exon3/6:c.322-3delT:,APTX:NM_001195250.1:splice_region&intron:LOW:exon3/6:c.364-3delT:,APTX:NM_001195251.1:splice_region&intron:LOW:exon5/8:c.484-3delT:,APTX:NM_001195252.1:splice_region&intron:LOW:exon4/7:c.310-3delT:,APTX:NM_175069.2:splice_region&intron:LOW:exon4/7:c.526-3delT:,APTX:NM_175073.2:splice_region&intron:LOW:exon5/8:c.484-3delT:
check(convert_hgvs2genomic("GRCh37", "NM_001195248.1", "c.526-3delT"),array("chr9",32986031,32986031,"A","-"));
#+-strand, splicing: chr17	18205751	18205760	GGAGAGTGAA	-	het	QUAL=781;DP=116;AF=0.32;MQM=60	TOP3A	splice_region&intron	TOP3A:NM_004618.3:splice_region&intron:LOW:exon6/18:c.644-12_644-3delTTCACTCTCC:
check(convert_hgvs2genomic("GRCh37", "NM_004618.3", "c.644-12_644-3delTTCACTCTCC"),array("chr17",18205751,18205760,"GGAGAGTGAA","-"));
#--strand, splicing: chr18	21111528	21111528	na	na	NPC1	c.3837+634_3837+637del	NM_000271.4
check(convert_hgvs2genomic("GRCh37", "NM_000271.4", "c.3837+634_3837+637del"),array("chr18",21111529,21111532,"CTTT","-"));
#--strand,splicing:	chr18	21112158	21112158	C	A	NPC1	c.3838+8G>T
check(convert_hgvs2genomic("GRCh37", "NM_000271.4", "c.3744_3747delCAGT"),array("chr18",21113326,21113329,"ACTG","-"));
#
check(convert_hgvs2genomic("GRCh37", "NM_004618.3","c.644-12del16"),array("chr17",18205745,18205760,"GCTCCTGGAGAGTGAA","-"));


//Dup
#chr1	54605318	54605318	-	G	het	QUAL=390;DP=18;AF=0.56;MQM=60	CDCP2	frameshift	CDCP2:NM_201546.3:frameshift:HIGH:exon4/4:c.1224dupC:p.Met409fs
check(convert_hgvs2genomic("GRCh37", "NM_201546.3", "c.1224dupC"),array("chr1",54605318,54605318,"-","G"));
#--strand:	chr8	48805816	48805816	-	G	hom	QUAL=12693;DP=451;AF=0.95;MQM=60	PRKDC	frameshift	PRKDC:NM_006904.6:frameshift:HIGH:exon31/86:c.3729dupC:p.Phe1244fs,PRKDC:NM_001081640.1:frameshift:HIGH:exon31/85:c.3729dupC:p.Phe1244fs
//@TODO error in cDNA position/conversion? Conversion fails... check(convert_hgvs2genomic("GRCh37", "NM_006904.6", "c.3729dupC"),array("chr8",48805816,48805816,"-","G"));
#+-strand: chr17	48452978	48452978	-	AGC	het	QUAL=9702;DP=943;AF=0.38;MQM=60	EME1	disruptive_inframe_insertion	EME1:NM_001166131.1:disruptive_inframe_insertion:MODERATE:exon2/9:c.410_412dupAGC:p.Lys137_Pro138insGln,EME1:NM_152463.2:disruptive_inframe_insertion:MODERATE:exon2/9:c.410_412dupAGC:p.Lys137_Pro138insGln
check(convert_hgvs2genomic("GRCh37", "NM_001166131.1", "c.410_412dupAGC"),array("chr17",48452978,48452978,"-","AGC"));
#WARNING: 'Start of variant in row 73 does not match converted start (c.3570_3573dupACTT, converted: chr18:21114432-21114432, previous: chr18:21114431-21114431). Skipping.' in /mnt/SRV017/users/ahschrc1/sandbox/NPC1_miriam/combine_annotate.php:94.
check(convert_hgvs2genomic("GRCh37", "NM_000271.4", "c.3570_3573dupACTT"),array("chr18",21114431,21114431,"-","AAGT"));


//@TODO NPC1:
#WARNING: 'Start of variant in row 144 does not match converted start (c.3245+1dupG, converted: chr18:21116638-21116637, previous: chr18:21115664-21115664). Skipping.' in /mnt/SRV017/users/ahschrc1/sandbox/NPC1_miriam/combine_annotate.php:94.
check(convert_hgvs2genomic("GRCh37", "NM_007294.3", "c.1A>G"),array("chr17",41276113,41276113,"T","C"));

end_test();

//##################################################################################
start_test("indel_for_vcf");
check(indel_for_vcf("GRCh37", "chr12",25368401,"TC","AG"),array("chr12",25368401,"TC","AG"));
check(indel_for_vcf("GRCh37", "chr12",25368401,"TC","-"),array("chr12",25368400,"TTC","T"));
check(indel_for_vcf("GRCh37", "chr12",25368401,"-","AG"),array("chr12",25368401,"T","TAG"));
end_test();

//##################################################################################
start_test("relative_path");
check(relative_path(src_folder(), data_folder()),"../test/data");
check(relative_path(data_folder(), src_folder()),"../../src");
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
check(vcf_varscan2("GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR","0/1:255:1560:1560:1070:489:31.37\%:2.7462E-165:36:35:707:363:355:134"), array(1560, 0.3137));
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

?>
