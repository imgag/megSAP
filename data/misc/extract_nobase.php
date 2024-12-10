<?php

include("../../src/Common/all.php");

$fasta = genome_fasta($argv[1]);

list($stdout) = execApptainer("ngs-bits", "FastaInfo", "-in $fasta -write_n nobase_regions.bed", [$fasta]);
print_r($stdout);

?>

