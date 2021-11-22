<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

/*
   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  463   1.3  0.6  1.7  chr1        10001   10468 (249240153) +  (TAACCC)n      Simple_repeat            1  463    (0)      1
 3612  11.4 21.5  1.3  chr1        10469   11447 (249239174) C  TAR1           Satellite/telo       (399) 1712    483      2
  484  25.1 13.2  0.0  chr1        11505   11675 (249238946) C  L1MC5a         LINE/L1             (2382) 5648   5452      3
  239  29.4  1.9  1.0  chr1        11678   11780 (249238841) C  MER5B          DNA/hAT-Charlie       (74)  104      1      4
  318  23.0  3.7  0.0  chr1        15265   15355 (249235266) C  MIR3           SINE/MIR             (119)  143     49      5
   18  23.2  0.0  2.0  chr1        15798   15849 (249234772) +  (TGCTCC)n      Simple_repeat            1   51    (0)      6
   18  13.7  0.0  0.0  chr1        16713   16744 (249233877) +  (TGG)n         Simple_repeat            1   32    (0)      7
  239  33.8 12.9  0.0  chr1        18907   19048 (249231573) +  L2a            LINE/L2               2942 3104  (322)      8
  994  31.2  6.0  2.5  chr1        19972   20405 (249230216) +  L3             LINE/CR1              2680 3129  (970)      9
  270  33.1  0.7  2.7  chr1        20531   20679 (249229942) +  Plat_L3        LINE/CR1              2802 2947  (639)     10
*/

$handle = fopen2("php://stdin", "r");
while (!feof($handle))
{
	//load
	$line = trim(fgets($handle));
	if ($line=="") continue;
	
	$line = preg_replace('/\s+/', ' ', $line);
	$parts = explode(" ", $line);
	if($parts[0]=="SW" || $parts[0]=="score") continue;
	
	$chr = $parts[4];
	$start = $parts[5]-1;
	$end = $parts[6];
	$repeat = $parts[9];
	$class = $parts[10];
	
	//write
	print "$chr\t$start\t$end\t{$repeat} ($class)\n";
}

fclose($handle);

?>