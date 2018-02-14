<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_Hgvs2Genomic", "Convert HGSV coordinates to genomic positions. Prints convertes position to stdout.");
$parser->addString("refseq",  "Refseq transcript ID.", false);
$parser->addString("hgvs",  "cDNA position in HGVS format.", false);
extract($parser->parse($argv));

list($chr,$start,$end,$ref,$obs) = convert_hgvs2genomic($refseq, $hgvs);
print "$chr\t$start\t$end\t$ref\t$obs";