<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("20260325_normalize_BNDs_in_NGSD", "Reorders the positions of translocations in NGSD so that pos1 <= pos2");
$parser->addString("ps", "Processed sample which should be updated.", false);
extract($parser->parse($argv));

$db = DB::getInstance("NGSD");

//get SV callset
$ps_info = get_processed_sample_info($db,$ps);
$ps_id = $info1['ps_id'];
$sv_callset_id = $db->getValue("SELECT id FROM sv_callset WHERE processed_sample_id = $ps_id");

//get SV ids which have to be updated
// callset id NA12878_45:	57136
// callset id DNA2405752A1:	51909 (contains BNDs on same chr)
$ids2update = $db->getValues("SELECT id FROM sv_translocation WHERE sv_callset_id={$sv_callset_id} AND ( (CAST(chr1 AS UNSIGNED) > CAST(chr2 AS UNSIGNED)) OR ((chr1 = chr2) AND (start1 > start2)) OR ((chr1 = chr2) AND (start1 = start2) AND (end1 > end2)) )");

//update SVs
$db->beginTransaction();
$hash_select = $db->prepare("SELECT chr1, start1, end1, chr2, start2, end2 FROM sv_translocation WHERE id=:id");
$hash_update = $db->prepare("UPDATE sv_transaction SET chr1=:chr1, start1=:start1, end1=:end1, chr2=:chr2, start2=:start2, end2=:end2 WHERE id=:id");
foreach ($ids2update as $sv_id) 
{
	//get values:
	$db->bind($hash_select, "id", $sv_id);
	$result = $db->execute($hash_select, true);
	
	$chr1 = $result[0]['chr1'];
	$start1 = $result[0]['start1'];
	$end1 = $result[0]['end1'];
	$chr2 = $result[0]['chr2'];
	$start2 = $result[0]['start2'];
	$end2 = $result[0]['end2'];

	//prepare update query (switched pos1/pos2)
	$db->bind($hash_update, "id", $sv_id);
	$db->bind($hash_update, "chr1", $chr2);
	$db->bind($hash_update, "start1", $start2);
	$db->bind($hash_update, "end1", $end2);
	$db->bind($hash_update, "chr2", $chr1);
	$db->bind($hash_update, "start2", $start1);
	$db->bind($hash_update, "end2", $end1);

	$db->execute($hash_update, true);
}

$db->endTransaction();
?>