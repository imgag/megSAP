<?php

/**
  @page vc_straglr

*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_straglr", "Call repeat expansions with straglr. Creates an BED file.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addOutfile("out", "Output BED file.", false);
$parser->addInfile("loci", "BED file containing repeat loci.", false);
//optional
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

// use BAM file name as fallback if no processed sample name is provided
if(!isset($pid)) $pid = basename($in);

//init
$out_prefix = dirname($out)."/".basename($out, ".vcf");
$straglr = get_path("straglr");


// prepare command
$args = [];
$args[] = $straglr;
$args[] = "--loci {$loci}";
$args[] = "--nprocs $threads";
$args[] = "$in";
$args[] = genome_fasta($build);
$args[] = $out_prefix;
		

// run straglr
$parser->exec(get_path("python3"), implode(" ", $args));


//TODO: Post-processing


?>