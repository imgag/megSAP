##fileformat=VCFv4.2
##INFO=<ID=REVEL,Number=1,Type=Float,Description="Maximum REVEL score of all scored transcripts.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$last_tag = "";
$max_score = -1;
$handle = fopen2("php://stdin", "r");
while($line = fgets($handle))
{
	$line = trim($line);
	if ($line=='' || $line[0]=='#') continue;
	
	list($chr, $pos_hg19, $pos, $ref, $alt, $aaref, $aaalt, $score) = explode(',', $line);
	if ($chr=='chr') continue;
	if ($pos=='.') continue;
	
	$tag = $chr.' '.$pos.' '.$ref.' '.$alt;
	if ($last_tag==$tag) //same variant
	{
		$max_score = max($max_score, $score);
	}
	else
	{
		//print out last variant
		if ($last_tag!='')
		{
			list($chr2, $pos2, $ref2, $alt2) = explode(' ', $last_tag);
			print "{$chr2}\t{$pos2}\t.\t{$ref2}\t{$alt2}\t30\t.\tREVEL={$max_score}\n";
		}
		
		//init for next line
		$last_tag = $tag;
		$max_score = $score;
	}
}

//print out last variant
list($chr2, $pos2, $ref2, $alt2) = explode(' ', $last_tag);
print "{$chr2}\t{$pos2}\t.\t{$ref2}\t{$alt2}\t30\t.\tREVEL={$max_score}\n";

//close handle
fclose($handle);

?>