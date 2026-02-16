<?php

/** 
	@page mapping_blatlike
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_blatlike", "Maps reads to a reference genome using bwa.");
$parser->addString("seq", "DNA sequence to map to the reference genome.", false);
//optional
$parser->addOutfile("out",  "Output TSV file. Writes to STDOUT if unset.", true);
$parser->addInt("secondary", "Maximum number of secondary alignments to keep.", true, 20);
$parser->addString("build", "Genome build to use.", true, "GRCh38");
$parser->addFlag("debug", "Enable debug output to stdout.");
extract($parser->parse($argv));

//init
$seq = strtoupper(trim($seq));
if (strlen($seq)<20) trigger_error("Sequence is shorter than 20bp. This is not supported!", E_USER_ERROR);
if (strlen($seq)>500) trigger_error("Sequence is longer than 500bp. This is not supported!", E_USER_ERROR);
if (!preg_match("/^[ACGT]+$/", $seq)) trigger_error("Sequence '$seq' contains not only A,C,G,T!", E_USER_ERROR);
$genome_fasta = genome_fasta($build);

//store sequence as FASTA
$tmp_fasta = $parser->tempFile(".fa");
file_put_contents($tmp_fasta, ">MY_READ\n{$seq}");

//run minimap
$args = [
	"mem",
	$genome_fasta,
	"-a",
	"-c {$secondary}",
	$tmp_fasta
];
$suffix = trim(get_path("bwa_mem2_suffix", false));
list($stdout, $stderr) = $parser->execApptainer("bwa-mem2", "bwa-mem2{$suffix}", implode(" ", $args), [$genome_fasta], []);

//debug output
if ($debug)
{
	foreach($stderr as $line)
	{
		print "##STRERR: ".trim($line)."\n";
	}
}

//parse hits
function match_length($cigar)
{
	$length = 0;
	
    preg_match_all('/(\d+)([MID=X])/i', $cigar, $matches_arr, PREG_SET_ORDER);
    foreach ($matches_arr as $m)
	{
		$op  = $m[2];
		if ($op=='M' || $op=='D') $length += $m[1];
    }

	return $length;
}

$hits = [];
foreach($stdout as $line)
{
	if (!starts_with($line, "MY_READ\t")) continue;
	
	if ($debug) print "READ: $line\n";
	
	$parts =  explode("\t", $line);
	list($id, $flags, $chr, $start, $map_q, $cigar) = $parts;
	
	//skip unmapped
	if ($chr=='*') continue;

	$length = match_length($cigar);
	$edit_distance = "";
	for($i=11; $i<count($parts); ++$i)
	{
		if (starts_with($parts[$i], "NM:i:")) $edit_distance = substr($parts[$i], 5);
	}
	$hits[] = [$chr.":".$start."-".($start+$length), $map_q, $length, $edit_distance];
}

//sort by mapping quality and match length
usort($hits, function ($a, $b)
	{
		return $b[1] <=> $a[1]  // mapping quality
			?: $b[2] <=> $a[2]; // length
	});

//make sure we only show the to X hits
$hits = array_slice($hits, 0, $secondary+1);

//write output
$h = fopen2($out=="" ? "php://stdout" : $out, "w");
fputs($h, "#position\tmapping_quality\tmatch_length\tedit_distance\n");
foreach($hits as $hit)
{
	fputs($h, implode("\t", $hit)."\n");
}
fclose($h);

?>