<?php 
/** 
	@page db_import_variants
	@todo fix somatic variant import (s. column indices)
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_variants", "Imports variants to NGSD.");
$parser->addString("id", "Processing ID (e.g. GS000123_01 for germline variants, GS000123_01-GS000124_01 for tumor-normal pairs).", false);
$parser->addInfile("var",  "Input variant list in TSV format.", false);
// optional
$parser->addEnum("mode",  "Import mode.", true, array("germline", "somatic"), "germline");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addFlag("force", "Overwrites already existing DB entries instead of throwing an error.");
$parser->addString("build", "The genome build to use.", true, "GRCh37");
extract($parser->parse($argv));

function getVariant($db, $id)
{
	$tmp = $db->executeQuery("SELECT chr, start, end, ref, obs FROM variant WHERE id='$id'");
	if (count($tmp)==0) trigger_error("No variant with id '$id' in database!", E_USER_ERROR);
	$var = array_values($tmp[0]);
	return $var[0].":".$var[1]."-".$var[2]." ".$var[3].">".$var[4];
}

//check pid format
$samples = array();
if(preg_match("/^([A-Za-z0-9]{4,})_(\d{2})$/", $id, $matches) && $mode="germline")
{
	$samples["germline"] = array("name"=>($matches[1]."_".$matches[2]), "sname"=>$matches[1],"pname"=>(int)$matches[2]);
}
else if(preg_match("/^([A-Za-z0-9]{4,})_(\d{2})-([A-Za-z0-9]{4,})_(\d{2})$/", $id, $matches) && $mode=="somatic")
{
	$samples["tumor"] = array("name"=>($matches[1]."_".$matches[2]),"sname"=>$matches[1],"pname"=>(int)$matches[2]);
	$samples["germline"] = array("name"=>($matches[3]."_".$matches[4]),"sname"=>$matches[3],"pname"=>(int)$matches[4]);
}
else trigger_error("'$id' is not a valid processing ID!", E_USER_ERROR);

//load variant file
$file = Matrix::fromTSV($var);

//database connection
$db_connect = DB::getInstance($db);

//convert GS-numbers to ids
$hash = $db_connect->prepare("SELECT s.id as sid, ps.id as psid FROM sample as s, processed_sample as ps WHERE s.id = ps.sample_id and s.name = :sname and ps.process_id = :pname;");
foreach($samples as $key => $sample)
{
	$db_connect->bind($hash, 'sname', $sample["sname"]);
	$db_connect->bind($hash, 'pname', $sample["pname"]);
	$db_connect->execute($hash, true); 
	$result = $db_connect->fetch($hash);
	if(count($result)!=1)
	{
		trigger_error("PID '".$samples[$key]["name"]."' not found in DB. Sample and Processed Sample entered to DB?", E_USER_ERROR);
	}
	$samples[$key]["sid"] = $result[0]['sid'];
	$samples[$key]["pid"] = $result[0]['psid'];
}
$db_connect->unsetStmt($hash);

// deal with germline variants
if($mode=="germline")
{
	$psid = $samples["germline"]["pid"];
	$sid = $samples["germline"]["sid"];
	
	//prevent tumor samples from being imported into the germline variant table
	if ($db_connect->getValue("SELECT tumor FROM sample WHERE id='$sid'")==1)
	{
		trigger_error("Skipping import into germline variant list because sample ".$samples["germline"]["name"]." is a tumor sample.", E_USER_WARNING);
		return;
	}
	
	//check if variants were already imported for this PID
	$count_old = $db_connect->getValue("SELECT count(*) FROM detected_variant WHERE processed_sample_id='".$psid."'");
	if($count_old!=0 && !$force)
	{
		trigger_error("Variants were already imported for PID '$id'. Use the flag '-force' to overwrite them.", E_USER_ERROR);
	}
	$parser->log("Found $count_old variants for PID '$id' already in DB.");

	$parser->log("DELETING DETECTED VARIANTS");
	$t_start = microtime(true);
	
	//remove old variants (and store class4/5 variants for check)
	$v_class = array(); //variant id => classification
	if ($count_old!=0 && $force)
	{
		//preserve classification
		$tmp = $db_connect->executeQuery("SELECT vc.variant_id, vc.class FROM detected_variant dv, variant_classification vc WHERE dv.processed_sample_id='$psid' AND dv.variant_id=vc.variant_id AND (vc.class='4' OR vc.class='5')");
		foreach($tmp as $row)
		{
			$v_class[$row['variant_id']] = $row['class'];
		}
		
		//remove old variants
		$db_connect->executeStmt("DELETE FROM detected_variant WHERE processed_sample_id='$psid'");
		$parser->log("Deleted old variants because the flag '-force' was used.");
	}

	//abort if there are no variants in the input file
	if ($file->rows()==0)
	{
		$parser->log("No variants imported (empty variant TSV file).");
		exit();
	}

	//get column indices of input file
	$i_chr = $file->getColumnIndex("chr");
	$i_sta = $file->getColumnIndex("start");
	$i_end = $file->getColumnIndex("end");
	$i_ref = $file->getColumnIndex("ref");
	$i_obs = $file->getColumnIndex("obs");
	$i_snp = $file->getColumnIndex("dbSNP");
	$i_10g = $file->getColumnIndex("1000g");
	$i_exa = $file->getColumnIndex("ExAC");
	$i_gno = $file->getColumnIndex("gnomAD");
	$i_gen = $file->getColumnIndex("gene");
	$i_typ = $file->getColumnIndex("variant_type");
	$i_cod = $file->getColumnIndex("coding_and_splicing");
	$i_geno = $file->getColumnIndex("genotype");
	
	$parser->log(" done ".time_readable(microtime(true) - $t_start));
	$parser->log("INSERT/UPDATE VARIANT");
	$c_ins = 0;
	$c_upd = 0;
	$t_start = microtime(true);
	
	//insert/update table 'variant'
	$var_ids = array();
	$hash = $db_connect->prepare("INSERT INTO variant (chr, start, end, ref, obs, dbsnp, 1000g, exac, gnomad, gene, variant_type, coding) VALUES (:chr, :start, :end, :ref, :obs, :dbsnp, :1000g, :exac, :gnomad, :gene, :variant_type, :coding) ON DUPLICATE KEY UPDATE id=id");
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);

		//remove GSvar additional information from string
		$dbsnp = trim($row[$i_snp]);
		if (contains($dbsnp, "[")) $dbsnp = substr($dbsnp, 0, strpos($dbsnp, "[")-1);
		
		//determine variant ID
		$variant_id = $db_connect->getValue("SELECT id FROM variant WHERE chr='".$row[$i_chr]."' AND start='".$row[$i_sta]."' AND end='".$row[$i_end]."' AND ref='".$row[$i_ref]."' AND obs='".$row[$i_obs]."'", -1);
		if ($variant_id!=-1) //update (common case)
		{
			//compose array of meta data
			$metadata = array(
							"dbsnp" => $dbsnp, 
							"1000g" => $row[$i_10g],
							"exac" => $row[$i_exa],
							"gnomad" => $row[$i_gno],
							"gene" => $row[$i_gen],
							"variant_type" => $row[$i_typ],
							"coding" => $row[$i_cod]
							);
							
			//check which variant meta data in NGSD needs to be updated
			$metadata_changed = array();
			$res = $db_connect->executeQuery("SELECT ".implode(", ", array_keys($metadata))." FROM variant WHERE id='".$variant_id."'");
			foreach($metadata as $key => $value)
			{
				if ($res[0][$key]!=$value)
				{
					$metadata_changed[$key] = $value;
				}
			}
			
			//update only the data that needs updating
			if (count($metadata_changed)!=0)
			{
				//prepare query
				$parts = array();
				foreach($metadata_changed as $key => $value)
				{
					$parts[] = "$key=:$key";
				}
				$hash2 = $db_connect->prepare("UPDATE variant SET ".implode(", ", $parts)." WHERE id=:id");
				
				//bind and execute
				foreach($metadata_changed as $key => $value)
				{
					$db_connect->bind($hash2, $key, $value, array(""));
				}
				$db_connect->bind($hash2, "id", $variant_id, array(""));
				$db_connect->execute($hash2, true);
				$db_connect->unsetStmt($hash2);
				++$c_upd;
			}
		}
		else //insert (rare case)
		{
			$db_connect->bind($hash, "chr", $row[$i_chr]);
			$db_connect->bind($hash, "start", $row[$i_sta]);
			$db_connect->bind($hash, "end", $row[$i_end]);
			$db_connect->bind($hash, "ref", $row[$i_ref]);
			$db_connect->bind($hash, "obs", $row[$i_obs]);
			$db_connect->bind($hash, "dbsnp", $dbsnp, array(""));
			$db_connect->bind($hash, "1000g", $row[$i_10g], array(""));
			$db_connect->bind($hash, "exac", $row[$i_exa], array(""));
			$db_connect->bind($hash, "gnomad", $row[$i_gno], array(""));
			$db_connect->bind($hash, "gene", $row[$i_gen]);
			$db_connect->bind($hash, "variant_type", $row[$i_typ], array(""));
			$db_connect->bind($hash, "coding", $row[$i_cod], array(""));
			$db_connect->execute($hash, true);
			++$c_ins;

			$variant_id = $db_connect->getValue("SELECT id FROM variant WHERE chr='".$row[$i_chr]."' AND start='".$row[$i_sta]."' AND end='".$row[$i_end]."' AND ref='".$row[$i_ref]."' AND obs='".$row[$i_obs]."'", -1);
		}
		$var_ids[] = $variant_id;
	}
	$db_connect->unsetStmt($hash);
	
	$parser->log(" inserted variants: $c_ins");
	$parser->log(" updated variants: $c_upd");
	$parser->log(" done ".time_readable(microtime(true) - $t_start));
	
	$parser->log("INSERT DETECTED VARIANT");
	$t_start = microtime(true);

	//insert variants into table 'detected_variant'
	$db_connect->beginTransaction();
	$hash = $db_connect->prepare("INSERT INTO detected_variant (processed_sample_id, variant_id, genotype) VALUES ('$psid', :variant_id, :genotype);");
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);
		$variant_id = $var_ids[$i];

		//skip invalid variants
		if ($row[$i_typ]=="invalid") continue;
		
		//remove class 4/5 variant from list (see check below)
		if (isset($v_class[$variant_id]))
		{
			unset($v_class[$variant_id]);
		}
		
		//bind
		$db_connect->bind($hash, "variant_id", $variant_id);
		$db_connect->bind($hash, "genotype", $row[$i_geno]);

		//execute
		$db_connect->execute($hash, true);
	}
	$db_connect->endTransaction();
	$db_connect->unsetStmt($hash);

	$parser->log(" done ".time_readable(microtime(true) - $t_start));
	
	//check that all important variant are still there (we unset all re-imported variants above)
	foreach($v_class as $v_id => $class)
	{
		trigger_error("Variant (".getVariant($db_connect, $v_id).") with classification '$class' is no longer in variant list!", E_USER_ERROR);
	}
}

if($mode=="somatic")
{
//somatic specific
	// check if variants were already imported for this PID
	$tmp = $db_connect->executeQuery("SELECT count(id) FROM detected_somatic_variant WHERE processed_sample_id_tumor='".$samples["tumor"]["pid"]."' AND processed_sample_id_normal='".$samples["germline"]["pid"]."'");
	$count_old = $tmp[0]['count(id)'];
	//throw error if not forced to overwrite and terms already exist
	if($count_old!=0 && !$force)
	{
		trigger_error("Variants were already imported for PIDs '".$samples["tumor"]["name"]."-".$samples["germline"]["name"]."'. Use the flag '-force' to overwrite them.", E_USER_ERROR);
	}
	$parser->log("Found $count_old variants for PID '".$samples["tumor"]["name"]."-".$samples["germline"]["name"]."' already in DB.");

	//remove old variants
	$comment_info = array();
	if ($count_old!=0 && $force)
	{
		//store validation info
		$comment_info = $db_connect->executeQuery("SELECT variant_id, comment FROM detected_somatic_variant WHERE processed_sample_id_tumor='".$samples["tumor"]["pid"]."' AND processed_sample_id_normal='".$samples["germline"]["pid"]."' AND comment IS NOT NULL");
		
		//remove old variants
		$db_connect->executeStmt("DELETE FROM detected_somatic_variant WHERE processed_sample_id_tumor='".$samples["tumor"]["pid"]."' AND processed_sample_id_normal='".$samples["germline"]["pid"]."'");
		$parser->log("Deleted old variants because the flag '-force' was used.");
	}
//somatic specific end

//no difference between germline and somatic
	//abort if there are no variants in the input file
	if ($file->rows()==0)
	{
		$parser->log("No variants imported (empty variant TSV file).");
		exit();
	}

	//get column indices of input file
	$i_chr = $file->getColumnIndex("chr");
	$i_sta = $file->getColumnIndex("start");
	$i_end = $file->getColumnIndex("end");
	$i_ref = $file->getColumnIndex("ref");
	$i_obs = $file->getColumnIndex("obs");
	$i_snp = $file->getColumnIndex("dbSNP");
	$i_10g = $file->getColumnIndex("1000g");
	$i_exa = $file->getColumnIndex("ExAC");
	$i_gno = $file->getColumnIndex("gnomAD");
	$i_gen = $file->getColumnIndex("gene");
	$i_typ = $file->getColumnIndex("variant_type");
	$i_cod = $file->getColumnIndex("coding_and_splicing");

	//get total number of variants before inserting
	$tmp = $db_connect->executeQuery("SELECT count(id) FROM variant;");
	$no_var_before = $tmp[0]['count(id)'];

	//insert variants into table 'variant'
	$hash = $db_connect->prepare("INSERT INTO variant (chr, start, end, ref, obs, dbsnp, 1000g, exac, gnomad, gene, variant_type, coding) VALUES (:chr, :start, :end, :ref, :obs, :dbsnp, :1000g, :exac, :gnomad, :gene, :variant_type, :coding)");
	$hash2 = $db_connect->prepare("UPDATE variant SET dbsnp=:dbsnp, 1000g=:1000g, exac=:exac, gnomad=:gnomad , gene=:gene, variant_type=:variant_type, coding=:coding WHERE id=:id");
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);

		//remove GSvar additional information from   string
		$dbsnp = trim($row[$i_snp]);
		if (contains($dbsnp, "[")) $dbsnp = substr($dbsnp, 0, strpos($dbsnp, "[")-1);
		
		//if variant already present => update details
		$variant_id = $db_connect->getValue("SELECT id FROM variant WHERE chr='".$row[$i_chr]."' AND start='".$row[$i_sta]."' AND end='".$row[$i_end]."' AND ref='".$row[$i_ref]."' AND obs='".$row[$i_obs]."'", -1);
		if ($variant_id!=-1)
		{
			$db_connect->bind($hash2, "dbsnp", $dbsnp, array(""));
			$db_connect->bind($hash2, "1000g", $row[$i_10g], array(""));
			$db_connect->bind($hash2, "exac", $row[$i_exa], array(""));
			$db_connect->bind($hash2, "gnomad", $row[$i_gno], array(""));
			$db_connect->bind($hash2, "gene", $row[$i_gen]);
			$db_connect->bind($hash2, "variant_type", $row[$i_typ], array(""));
			$db_connect->bind($hash2, "coding", $row[$i_cod], array(""));
			$db_connect->bind($hash2, "id", $variant_id);
			$res = $db_connect->execute($hash2, true);
			continue;
		}

		//insert new variant		
		$db_connect->bind($hash, "chr", $row[$i_chr]);
		$db_connect->bind($hash, "start", $row[$i_sta]);
		$db_connect->bind($hash, "end", $row[$i_end]);
		$db_connect->bind($hash, "ref", $row[$i_ref]);
		$db_connect->bind($hash, "obs", $row[$i_obs]);
		$db_connect->bind($hash, "dbsnp", $dbsnp, array(""));
		$db_connect->bind($hash, "1000g", $row[$i_10g], array(""));
		$db_connect->bind($hash, "exac", $row[$i_exa], array(""));
		$db_connect->bind($hash, "gnomad", $row[$i_gno], array(""));
		$db_connect->bind($hash, "gene", $row[$i_gen]);
		$db_connect->bind($hash, "variant_type", $row[$i_typ], array(""));
		$db_connect->bind($hash, "coding", $row[$i_cod], array(""));

		//execute
		$db_connect->execute($hash, true);
	}

	//log the number of new variants
	$tmp = $db_connect->executeQuery("SELECT count(id) FROM variant;");
	$no_new_variants = $tmp[0]['count(id)'] - $no_var_before;
	$parser->log("Inserted $no_new_variants rows into table 'variant'.");
//general end

//somatic specific
	//insert variants into table 'detected_somatic_variant'
	$rescued_val_count = 0;
	$rescued_comment_count = 0;
	$hash = $db_connect->prepare("INSERT INTO detected_somatic_variant (processed_sample_id_tumor,processed_sample_id_normal,variant_id,variant_frequency,depth,quality_snp,comment) VALUES ('".$samples["tumor"]["pid"]."','".$samples["germline"]["pid"]."',:variant_id,:variant_frequency,:depth,:snp_q,:comment);");
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);

		//get variant_id
		$variant_id = $db_connect->getValue("SELECT id FROM variant WHERE chr='".$row[$i_chr]."' AND start='".$row[$i_sta]."' AND end='".$row[$i_end]."' AND ref='".$row[$i_ref]."' AND obs='".$row[$i_obs]."'", -1);

		//get comment from previous analysis (if it was available)
		$comment = NULL;
		for($v=0; $v<count($comment_info); ++$v)
		{
			$entry = $comment_info[$v];
			if ($entry["variant_id"]==$variant_id)
			{
				$comment = $entry["comment"];
				++$rescued_comment_count;
				array_splice($comment_info, $v, 1); //remove to check how many comment are missing after re-import (see below)
				break;
			}
		}

		//bind
		$i_frq = $file->getColumnIndex("tumor_af");
		$i_dep = $file->getColumnIndex("tumor_dp");
		$i_qsn = $file->getColumnIndex("snp_q");
		$db_connect->bind($hash, "variant_id", $variant_id);
		$db_connect->bind($hash, "variant_frequency", $row[$i_frq], array("n/a"));
		$db_connect->bind($hash, "depth", $row[$i_dep]);
		$db_connect->bind($hash, "snp_q", $row[$i_qsn], array("."));
		$db_connect->bind($hash, "comment", $comment);
		$db_connect->execute($hash, true);

//end somatic specific
	}
}


?>
