<?php

include("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("gsvar2bed", "Determines regions with high error rates.");
$parser->addInfileArray("in", "Input GSvar files.", false);
$parser->addOutfile("out", "Output BED file.", false);
//optional
$parser->addInt("errors", "Minimum number of samples that have errors in a region.", true, 2);
extract($parser->parse($argv));

//process input files
$c_gsvar = 0;
$regs = [];
foreach($in as $gsvar)
{
	$gsvar = trim($gsvar);
	if ($gsvar=="") continue;
	
	//print "processing: ".basename($gsvar)."\n";
	++$c_gsvar;
	$file = file($gsvar);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		list($chr, $start, $end, $ref, $obs) = explode("\t", $line);
		$var = "$chr:$start-$end $ref>$obs";
		if ($ref=="-") $end += 1; //insersion (two bases around insertion)
		if ($obs=="-") $start -= 1; //deletion (base before deletion)
		for($p=$start; $p<=$end; ++$p)
		{
			if (!isset($regs[$chr][$p]))
			{
				$regs[$chr][$p] = 0;
			}
			$regs[$chr][$p] += 1;
		}
	}
}

//write output BED file
$output = array();
foreach($regs as $chr => $tmp)
{
	foreach($tmp as $pos => $count)
	{
		if ($count>=$errors)
		{
			$output[] = "$chr\t".($pos-1)."\t".($pos)."\t$count\n";
		}
	}
}
file_put_contents($out, $output);

//merge regions
$parser->exec(get_path("ngs-bits")."BedMerge", "-in $out -out $out", false);

//annotate with gene names
$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $out -clear -extend 20 -out $out", false);

//statistics (size)
list($stdout) = $parser->exec(get_path("ngs-bits")."BedInfo", "-in $out", false);
print $stdout[0]."\n";
print $stdout[1]."\n";

$regs_per_gene = [];
$file = file($out);
foreach($file as $line)
{
	$line = nl_trim($line);
	if ($line=="" || $line[0]=="#") continue;
	list($chr, $start, $end, $genes) = explode("\t", $line);
	$genes = explode(", ", $genes);
	foreach($genes as $gene)
	{
		if (!isset($regs_per_gene[$gene])) $regs_per_gene[$gene] = 0;
		$regs_per_gene[$gene] += 1;
	}
}
arsort($regs_per_gene, SORT_NUMERIC);
foreach($regs_per_gene as $gene => $count)
{
	print $gene."\t".$count."\n";
}
	
		


?>
