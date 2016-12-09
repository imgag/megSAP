<?php	

/** 
	@page plink_regions
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("plink_regions", "\$Rev: 712 $", "Extracts homocygosity regions from the PLINK output.");
$parser->addInfile("in",  "PLINK homocygosity mapping output file.", false);
$parser->addOutfile("out",  "Output BED file.", false);
$parser->addInt("gap",  "Gaps between regions upto this size (in bases) are automatically closed.", true, 1000000);
extract($parser->parse($argv));

// (1) calculate interesting regions (affected persons are homocygote, not-affected are not homocygote)
$regs_a = array();
$regs_na = array();
$file = Matrix::fromTSV($in);
for($i=0; $i<$file->rows(); ++$i)
{
	list(, $iid, $phe, $chr, , , $start, $end) = $file->getRow($i);
	
	if ($phe==2) //affected (hash by IID)
	{
		if (!isset($regs_a[$iid]))
		{
			$regs_a[$iid] = array();
		}
		$regs_a[$iid][] = "chr$chr\t$start\t$end";
	}
	else if ($phe==1) // not affected
	{
		$regs_na[] = "chr$chr\t$start\t$end";
	}
	else
	{
		trigger_error("Unknown phenotype '$phe' encountered!", E_USER_ERROR);
	}
}

// intersect homocygote regions of affected persons
$tmp1 = $parser->tempFile(".bed");
$tmp2 = $parser->tempFile(".bed");
$is_first = true;
foreach($regs_a as $iid => $regs)
{
	if ($is_first)
	{
		file_put_contents($tmp1, implode("\n", $regs));
		$is_first = false;
		continue;
	}
	
	file_put_contents($tmp2, implode("\n", $regs));
	$parser->exec(get_path("ngs-bits")."BedIntersect", "-in $tmp1 -in2 $tmp2 -out $tmp1", true);
}

// subtract homoygote regions of not affected persons
file_put_contents($tmp2, implode("\n", $regs_na));
$parser->exec(get_path("ngs-bits")."BedSubtract", "-in $tmp1 -in2 $tmp2 -out $tmp1", true);

// close gaps
$parser->exec(get_path("ngs-bits")."BedExtend", "-in $tmp1 -out $tmp1 -n ".($gap/2), true);
$parser->exec(get_path("ngs-bits")."BedMerge", "-in $tmp1 -out $tmp1", true);
$parser->exec(get_path("ngs-bits")."BedShrink", "-in $tmp1 -out $out -n ".($gap/2), true);

?>