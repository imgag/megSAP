<?php

include("../../src/Common/all.php");

$fasta = genome_fasta($argv[1]);

list($stdout) = exec2(get_path("ngs-bits")."FastaInfo -in $fasta -write_n nobase_regions.bed");
print_r($stdout);

?>

