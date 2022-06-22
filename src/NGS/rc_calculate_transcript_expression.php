<?php
/**
 * @page rc_calculate_transcript_expression
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parameters
$parser = new ToolBase("rc_calculate_transcript_expression", "Calculates transcript expression based on expression of exons.");
$parser->addInfile("in", "Input normalized exon read counts file.", false, true);
$parser->addOutfile("out", "Output file for normalized transcript read counts.", false);
$parser->addFlag("only_gencode_basic", "Use only Gencode Basic transcripts");
$parser->addFlag("only_ensembl_canonical", "Use only ensembl canonical transcripts");
$parser->addFlag("only_mane_select", "Use only Mane Select transcripts");
$parser->addFlag("only_mane_plus_clinical", "Use only Mane plus Clinical transcripts");


//optional parameters
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");

extract($parser->parse($argv));



//parse exon expression
$exons = array();

$fp = fopen2($in, "r");
while (($line = fgets($fp)) !== false) 
{
	//skip comment and header
	if($line[0] == "#") continue;
	
	$split_line = explode("\t", $line);
	$tmp = explode(":", $split_line[2]);
	$pos = $tmp[0].":".$tmp[1];
	if (!starts_with($pos, "chr")) $pos = "chr".$pos;
	$raw = intval($split_line[3]);
	$rpb = floatval($split_line[4]);
	$srpb = floatval($split_line[5]);

	$exons[$pos] = array($raw, $rpb, $srpb);
}

fclose($fp);


//get all transcripts

//establish database connection
$db = DB::getInstance($db);

//limit transcript set
$additional_conditions = array();
if($only_gencode_basic) $additional_conditions[] = "is_gencode_basic=1";
if($only_ensembl_canonical) $additional_conditions[] = "is_ensembl_canonical=1";
if($only_mane_select) $additional_conditions[] = "is_mane_select=1";
if($only_mane_plus_clinical) $additional_conditions[] = "is_mane_plus_clinical=1";

$additional_conditions_str = "";
if(count($additional_conditions) > 0) $additional_conditions_str = " WHERE ".implode(" AND ", $additional_conditions)." ORDER BY name";

//prepare queries
$query_transcript_info = "SELECT t.name, g.ensembl_id, g.symbol, t.biotype, t.chromosome FROM gene_transcript t INNER JOIN gene g ON g.id=t.gene_id WHERE t.id=:0";
$db->prepare($query_transcript_info);
$query_exons = "SELECT `start`, `end` FROM gene_exon WHERE transcript_id=:0";
$db->prepare($query_exons);
$query_transcripts = "SELECT id FROM gene_transcript".$additional_conditions_str;
$transcript_ids = $db->getValues($query_transcripts);


$output_rows = array();
$debug_exons = array();
$skipped_transcripts = 0;
$skipped_transcripts_ccds = 0;
$used_transcripts = 0;
$skipped_exons = 0;
$used_exons = 0;
$failed_transcripts = array();
foreach ($transcript_ids as $id) 
{
	//prepare vars
	$values_raw = array();
	$values_rpb = array();
	$values_srpb = array();

	//get transcript info
	$res = $db->executeQuery($query_transcript_info, ["0"=>$id]);
	$transcript_name = $res[0]["name"];
	$ensg = $res[0]["ensembl_id"];
	$symbol = $res[0]["symbol"];
	$biotype = $res[0]["biotype"];
	$chr = $res[0]["chromosome"];

	if(starts_with($transcript_name, "CCDS"))
	{
		$skipped_transcripts_ccds++;
		continue;
	}

	//get exons
	$res = $db->executeQuery($query_exons, ["0"=>$id]);
	$all_exons_found = true;
	foreach($res as $row)
	{
		$pos = "chr".$chr.":".$row["start"]."-".$row["end"];

		$debug_exons[] = $ensg."\t".$pos;

		//get expression data for exon
		if (!array_key_exists($pos, $exons))
		{
			user_error("WARNING: Exon $pos not found in exon file", E_USER_WARNING);
			$all_exons_found = false;
			$skipped_exons++;
			continue;
		}
		else
		{
			$values_raw[]  = $exons[$pos][0];
			$values_rpb[]  = $exons[$pos][1];
			$values_srpb[] = $exons[$pos][2];
			$used_exons++;
		}
	}

	//calculate combined stats
	$raw = array_sum($values_raw);
	$rpb = 0;
	if (count($values_rpb) > 0) $rpb = array_sum($values_rpb)/count($values_rpb);
	$srpb = 0;
	if(count($values_srpb) > 0) $srpb = array_sum($values_srpb)/count($values_srpb);
	if($all_exons_found)
	{
		$output_rows[] = array($transcript_name, $raw, $rpb, $srpb, $ensg, $symbol, $biotype);
		$used_transcripts++;
	}
	else
	{
		$failed_transcripts[] = $transcript_name;
		$skipped_transcripts++;
	}
	
}


//create output file
$fp = fopen2($out, "w");

//write comments
fwrite($fp, "##");

fwrite($fp, "##DESCRIPTION=transcript_id=Transcript transcript_name identifier\n");
fwrite($fp, "##DESCRIPTION=raw=Number of reads overlapping transcript exons\n");
fwrite($fp, "##DESCRIPTION=rpb=Reads overlapping transcript exon per base\n");
fwrite($fp, "##DESCRIPTION=srpb=Reads overlapping transcript exon per base, scaled by total number of reads in 100 million\n");
fwrite($fp, "##DESCRIPTION=gene_id=Gene ENSG identifier\n");
fwrite($fp, "##DESCRIPTION=symbol=Gene symbol\n");
fwrite($fp, "##DESCRIPTION=biotype=transcript biotype\n");

//write header
fwrite($fp, "#transcript_id\traw\trpb\tsrpb\tgene_id\tsymbol\tbiotype\n");

//write data
foreach ($output_rows as $row) 
{
	fwrite($fp, implode("\t", $row)."\n");
}

fclose($fp);

print "Exons used:    \t{$used_exons}\n";
print "Exons skipped: \t{$skipped_exons}\n";
print "Transcripts used:    \t{$used_transcripts}\n";
print "Transcripts skipped: \t{$skipped_transcripts}\n";
print "Transcripts skipped (CCDS): \t{$skipped_transcripts_ccds}\n";

#print "Skipped transcripts:\n".implode(", ", $failed_transcripts);

file_put_contents("db_exons.tsv", implode("\n", $debug_exons));

?>