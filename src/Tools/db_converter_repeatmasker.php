<?php

/*
##gff-version 3
##sequence-region chr1 1 249250621
chr1	RepeatMasker	dispersed_repeat	10001	10468	1504	+	.	Target=(CCCTAA)n 1 463
chr1	RepeatMasker	dispersed_repeat	10469	11447	3612	-	.	Target=TAR1 483 1712
chr1	RepeatMasker	dispersed_repeat	11505	11675	484	-	.	Target=L1MC5a 199 395
chr1	RepeatMasker	dispersed_repeat	11678	11780	239	-	.	Target=MER5B 1 104
chr1	RepeatMasker	dispersed_repeat	15265	15355	318	-	.	Target=MIR3 49 143
chr1	RepeatMasker	dispersed_repeat	16713	16749	203	+	.	Target=(TGG)n 1 37
chr1	RepeatMasker	dispersed_repeat	18907	19048	239	+	.	Target=L2a 2942 3104
chr1	RepeatMasker	dispersed_repeat	19972	20405	994	+	.	Target=L3 2680 3129
*/

$handle = fopen("php://stdin", "r");
while (!feof($handle))
{
	//load
	$line = trim(fgets($handle));
	if ($line=="" || $line[0]=="#") continue;
	
	list($chr, , ,$start, $end, , , , $details) = explode("\t", $line);
	
	//remove 'Target=' and numbers
	$details = substr($details, 7);
	$details = explode(" ", $details);
	$details = $details[0];
	
	//convert position (1-based, end included => 0-based, end excluded)
	$start -= 1;
	
	//write
	print "$chr\t$start\t$end\t$details\n";
}

fclose($handle);

?>