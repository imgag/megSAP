<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//create map: transcript DB ID > ENST 
//>zgrep 803086 transcript.txt.gz
//803086  156171  8355    131550  944203  959256  -1      3531052 ensembl_havana  protein_coding  \N      1       1306389 ENST00000327044 7       2008-04-29 11:17:41  2018-11-21 17:23:49
$transcript_id2enst = [];
$handle = gzopen2($argv[3], "r");
while (!feof($handle))
{
	//load
	$line = nl_trim(fgets($handle));
	if ($line=="") continue;

	$parts = explode("\t", $line);
	if (count($parts)!=17)
	{
		print_r($parts);
		exit(1);
	}
	
	$enst = $parts[13];
	if (!starts_with($enst, "ENST")) continue; //only Ensembl transcripts
	
	$transcript_id = $parts[0];
	
	if (isset($transcript_id2enst[$transcript_id])) trigger_error("transcript id to ENST mapping not unique!", E_USER_ERROR);
	$transcript_id2enst[$transcript_id] = $enst;
}
fclose($handle);

//create map transcript DB ID > translation DB ID
//>zgrep 1306389 translation.txt.gz
//1306389 803086  17      2172846 107     2172864 ENSP00000317992 6       2008-04-29 11:17:41     2009-08-05 14:27:16
$translation_id2transcript_id = [];
$handle = gzopen2($argv[2], "r");
while (!feof($handle))
{
	//load
	$line = nl_trim(fgets($handle));
	if ($line=="") continue;

	$parts = explode("\t", $line);
	if (count($parts)!=10)
	{
		print_r($parts);
		exit(1);
	}
	$translation_id = $parts[0];
	$transcript_id = $parts[1];
	if (isset($translation_id2transcript_id[$translation_id])) trigger_error("translation id to transcript id mapping not unique!", E_USER_ERROR);
	$translation_id2transcript_id[$translation_id] = $transcript_id;
}
fclose($handle);

//parse domains
//>zgrep PF03715 protein_feature.txt.gz
//protein_feature.txt.gz:46686559 1306389 328     622     2       298     PF03715 10258   397.7   9.2e-116        \N      \N      Noc2    \N      \N

print "#ENST\tAA_START\tAA_END\tPFAM_ID\tPFAM_DESC\n";
$c_skipped_no_transcript = 0;
$c_skipped_no_enst = 0;
$domains = [];
$handle = gzopen2($argv[1], "r");
while (!feof($handle))
{
	//load
	$line = nl_trim(fgets($handle));
	if ($line=="") continue;

	$parts = explode("\t", $line);
	if (count($parts)!=15)
	{
		print_r($parts);
		exit(1);
	}

	$pfam_id = $parts[6];
	if (!starts_with($pfam_id, "PF")) continue; //only PFAM domains
	$pf_desc = $parts[12];
	$aa_start = $parts[2];
	$aa_end = $parts[3];
	$translation_id = $parts[1];
	if (!isset($translation_id2transcript_id[$translation_id]))
	{
		++$c_skipped_no_transcript;
		continue;
	}
	$transcript_id = $translation_id2transcript_id[$translation_id];
	if (!isset($transcript_id2enst[$transcript_id]))
	{
		++$c_skipped_no_enst;
		continue;
	}
	$enst = $transcript_id2enst[$transcript_id];
	
	print implode("\t", [$enst, $aa_start, $aa_end, $pfam_id, $pf_desc])."\n";
}
fclose($handle);

print "##Skipped domains, because of no transcript: {$c_skipped_no_transcript}\n";
print "##Skipped domains, because of no transcript: {$c_skipped_no_enst}\n";

?>