<?php
/** 
	@page genome_statistics
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("genome_statistics", "Writes statistics for a genome in FASTA format.");
$parser->addInfile("in", "Input VCF file.", false);
extract($parser->parse($argv));

$bases = array("A", "C", "G", "T", "N");

//parse input
$handle = fopen2($in, 'r');
$stats = array();
while(!feof($handle))
{
	$line = trim(fgets($handle));
	if ($line=="") continue;
	
	if ($line[0]==">")
	{
		list($chr) = explode(" ", $line);
		$chr = strtr($chr, array(">"=>"", "chr"=>""));
		print "##Processing chromosome $chr...\n";
		
		$base_idx = 0;
		foreach($bases as $base)
		{
			$stats[$chr][$base] = 0;
		}
	}
	else
	{
		$line = strtoupper($line);
		for ($i=0; $i<strlen($line); ++$i)
		{
			if (!in_array($line[$i], $bases))
			{
				print "## - invalid base ".$line[$i]." (".$chr.":".$base_idx.")\n";
				$line[$i] = "N";
			}
			$stats[$chr][$line[$i]] += 1;
			++$base_idx;
		}
	}
}

//print output
print "#chr\t".implode("\t", $bases)."\tsize\n";
foreach($stats as $chr => $data)
{
	print $chr;
	foreach($bases as $base)
	{
		print "\t".$data[$base];
	}
	print "\t".array_sum(array_values($data))."\n";
}



?>

