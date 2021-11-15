<?php

include("framework.php");

//##################################################################################
start_test("get_genome_build");
check(get_genome_build(data_folder()."/get_genome_build_bwaGRCh37.bam"), "GRCh37");
check(get_genome_build(data_folder()."/get_genome_build_bwaGRCh38.bam"), "GRCh38");
check(get_genome_build(data_folder()."/get_genome_build_bwamem2GRCh37.bam"), "GRCh37");
check(get_genome_build(data_folder()."/get_genome_build_bwamem2GRCh38.bam"), "GRCh38");
check(get_genome_build(data_folder()."/get_genome_build_dragenGRCh37.bam"), "GRCh37");
check(get_genome_build(data_folder()."/get_genome_build_dragenGRCh38.bam"), "GRCh38");
end_test();


?>
