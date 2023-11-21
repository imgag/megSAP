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
$out_prefix = dirname($out)."/".basename($out, ".bed");
$straglr = get_path("straglr");


// prepare command
$args = [];
$args[] = $straglr;
$args[] = "--loci {$loci}";
$args[] = "--nprocs $threads";
$args[] = "$in";
$args[] = genome_fasta($build);
$args[] = $out_prefix;
		
// prepare PATH
putenv("PATH=".getenv("PATH").PATH_SEPARATOR.dirname(get_path("trf")).PATH_SEPARATOR.dirname(get_path("blastn")));
print "\n........\n".getenv("PATH")."\n........\n";
// run straglr
$parser->exec(get_path("python3"), implode(" ", $args));

// annotate repeat names
//read and index catalog
$loci_content = file($loci, FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
$catalog = array();
foreach ($loci_content as $line) 
{
  if (starts_with($line, "#")) continue;
  list($chr, $start, $end, $motive, $repeat_id, $associated_disorder, $ref_size, $normal_range, $pathogenic_range) = explode("\t", $line);
  $catalog["{$chr}:{$start}-{$end}"] = array($repeat_id, $associated_disorder, $ref_size, $normal_range, $pathogenic_range);
}
//annotate catalog to output file
$bed_content_in = file($out, FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
$bed_content_out = array();
foreach ($bed_content_in as $line) 
{
  if (starts_with($line, "##")) $bed_content_out[] = $line;
  elseif (starts_with($line, "#")) $bed_content_out[] = $line."\trepeat_id";
  else 
  {
    list($chr, $start, $end, $motive) = explode("\t", $line);
    $bed_content_out[] = $line."\t".implode("\t", $catalog["{$chr}:{$start}-{$end}"]);
  }
}
file_put_contents($out, implode("\n", $bed_content_out));


//TODO: further post-processing?
?>