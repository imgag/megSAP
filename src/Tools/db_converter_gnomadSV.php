<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$handle = fopen2("php://stdin", "r");
$header_parsed = false;
while($line = fgets($handle))
{
	$line = trim($line);
	if ($line=="") continue;
	if($line[0]=="#")
	{
		// skip comments
		if($line[1]=="#") continue;


		// parse header line
		$header_columns = explode("\t", substr($line, 1));

		$idx_chr1 = array_search("chrom", $header_columns);
		$idx_start1 = array_search("start", $header_columns);
		$idx_end1 = array_search("end", $header_columns);
		$idx_chr2 = array_search("CHR2", $header_columns);
		$idx_start2 = array_search("POS2", $header_columns);
		$idx_end2 = array_search("END2", $header_columns);
		$idx_svtype = array_search("SVTYPE", $header_columns);
		$idx_algorithm = array_search("ALGORITHMS", $header_columns);

		$idx_an = array_search("AN", $header_columns);
		$idx_af = array_search("AF", $header_columns);
		$idx_nhomalt = array_search("N_HOMALT", $header_columns);
		$idx_maleac = array_search("MALE_AC", $header_columns);
		$idx_par = array_search("PAR", $header_columns);


		// print comment


		// print header
		print "#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tTYPE\tALGORITHMS\tAN\tAF\tHEMI\tHOM\n";
		$header_parsed = true;
		continue;
	}

	if (!$header_parsed) trigger_error("No header found in file!", E_USER_ERROR);
	
	// parse content line
	$parts = explode("\t", $line);

	$chr1 = $parts[$idx_chr1];
	$start1 = $parts[$idx_start1];
	$end1 = $parts[$idx_end1];
	$chr2 = $parts[$idx_chr2];
	$start2 = $parts[$idx_start2];
	$end2 = $parts[$idx_end2];
	$svtype = $parts[$idx_svtype];
	$algorithm = $parts[$idx_algorithm];

	$an = $parts[$idx_an];
	$af = $parts[$idx_af];
	$nhomalt = $parts[$idx_nhomalt];
	$maleac = $parts[$idx_maleac];
	$par = $parts[$idx_par] == "True";

	// skip MCNVs 
	if ($svtype == "MCNV") continue;

	//get required infos
	$is_chrx = ($chr1=="X") || ($chr2=="X");

	if ($is_chrx && !$par) 
	{
		$hemi = $maleac;
	} 
	else
	{
		$hemi=".";
	} 

	$af = ($an<200 ? 0.0 : number_format($af, 5));
	
	// write to output
	print "$chr1\t$start1\t$end1\t$chr2\t$start2\t$end2\t$svtype\t$algorithm\t$an\t$af\t$hemi\t$nhomalt\n";
}
fclose($handle);

?>