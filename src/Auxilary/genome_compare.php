<?php
/** 
	@page genome_compare
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("genome_compare", "Compares genome FASTA files.");
$parser->addInfile("in1", "Input genome FASTA file.", false);
$parser->addInfile("in2", "Input genome FASTA file.", false);
extract($parser->parse($argv));

//compare chromosomes
print "### comparing chromosomes ###\n";
$chrs1 = genome_chr_sizes($in1);
$chrs2 = genome_chr_sizes($in2);
$chrs_common = [];
foreach($chrs1 as $chr => $size)
{
	if (!isset($chrs2[$chr]))
	{
		print "##Notice: $chr missing in $in2\n";
	}
	else
	{
		$chrs_common[$chr] = $size;
	}
}
foreach($chrs2 as $chr => $size)
{
	if (!isset($chrs1[$chr]))
	{
		print "##Notice: $chr missing in $in1\n";
	}
}

print "##Common chromosomes: ".implode(", ", array_keys($chrs_common))."\n";


//check chromosome sizes match
print "### comparing chromosome sizes ###\n";
foreach($chrs_common as $chr => $size)
{
	if ($chrs1[$chr]!=$chrs2[$chr])
	{
		trigger_error("Chromosome size mismatch for $chr: ".$chrs1[$chr]." vs ".$chrs2[$chr], E_USER_ERROR);
	}
}

//find differing bases
print "### comparing chromosome sequences ###\n";
$diffs = [];
foreach($chrs_common as $chr => $size)
{
	print "##diffing $chr...\n";
	list($output1) = $parser->execApptainer("samtools", "samtools faidx", "{$in1} {$chr}:1-{$size} 2>&1", [$in1]);
	$output1 = strtoupper(trim(implode("", array_slice($output1, 1))));
	list($output2) = $parser->execApptainer("samtools", "samtools faidx", "{$in2} {$chr}:1-{$size} 2>&1", [$in2]);
	$output2 = strtoupper(trim(implode("", array_slice($output2, 1))));
	if (strlen($output1)!=strlen($output2)) trigger_error("Differing sequence length in $chr: ".strlen($output1)." vs ".strlen($output2), E_USER_ERROR);
	if($output1!=$output2)
	{
		$diff_count_chr = 0;
		for($i=0; $i<strlen($output1); ++$i)
		{
			if ($output1[$i]!=$output2[$i])
			{
				$diffs[] = [$chr, $i, $output1[$i], $output2[$i]];
				++$diff_count_chr;
			}
		}
		if ($diff_count_chr>0) print "##Base differences on $chr: $diff_count_chr\n";
	}
}

//output
print "#chr\tstart\tend\tsize\t$in1\t$in2\n";
$current_chr = null;
$start = null;
$end = null;
$seq1 = "";
$seq2 = "";
foreach ($diffs as list($chr, $pos, $base1, $base2))
{
    //continue current block if position is directly adjacent
    if ($current_chr === $chr && $pos === $end + 1)
    {
        $end = $pos;
        $seq1 .= $base1;
        $seq2 .= $base2;
    }
    else
    {
        //print previous block
        if ($current_chr !== null)
        {
            print "$current_chr\t$start\t$end\t".($end-$start)."\t$seq1\t$seq2\n";
        }

        //start new block
        $current_chr = $chr;
        $start = $pos;
        $end = $pos;
        $seq1 = $base1;
        $seq2 = $base2;
    }
}

//print final block
if ($current_chr !== null)
{
    print "$current_chr\t$start\t$end\t".($end-$start)."\t$seq1\t$seq2\n";
}

?>
