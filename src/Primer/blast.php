<?php

/** 
	@page blast
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("blast", "Blast a FASTA file against the hg19 database.");
$parser->addInfile("in",  "Input FASTA file.", false);
$parser->addOutfile("out",  "Output BLAST result file.", false);
extract($parser->parse($argv));

list($output) = $parser->exec("blastn", "-version", false);
print "BLAST info: ".$output[0];

$parser->exec("blastn", "-task blastn-short -db ".get_path("GRCh37_data_folder")."/dbs/blast/hg19 -query $in -outfmt 7 -out $out -num_threads 8", true);


?>