<?php 
/** 
	@page db_import_variants
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_variants", "Imports variants to NGSD.");
$parser->addString("id", "Processing ID (e.g. GS000123_01 for germline variants, GS000123_01-GS000124_01 for tumor-normal pairs).", false);
$parser->addInfile("var",  "Input variant list in GSvar format.", false);
// optional
$parser->addEnum("mode",  "Import mode.", true, array("germline", "somatic"), "germline");
$parser->addFloat("max_af", "Maximum allele frequency for import of germline variants.", true, 0.05);
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addFlag("force", "Overwrites already existing DB entries instead of throwing an error.");
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//returns a human-readable variant representation, e.g. chr5:749999393 A>G
function getVariant($db, $id)
{
	$tmp = $db->executeQuery("SELECT chr, start, end, ref, obs FROM variant WHERE id='$id'");
	if (count($tmp)==0) trigger_error("No variant with id '$id' in database!", E_USER_ERROR);
	$var = array_values($tmp[0]);
	return $var[0].":".$var[1]."-".$var[2]." ".$var[3].">".$var[4];
}

//adds variants or updates variants
function updateVariantTable($parser, $db_connect, $file, $max_af=-1.0)
{	
	//get column indices of input file
	$i_chr = $file->getColumnIndex("chr");
	$i_sta = $file->getColumnIndex("start");
	$i_end = $file->getColumnIndex("end");
	$i_ref = $file->getColumnIndex("ref");
	$i_obs = $file->getColumnIndex("obs");
	$i_10g = $file->getColumnIndex("1000g");
	$i_gno = $file->getColumnIndex("gnomAD");
	$i_gen = $file->getColumnIndex("gene");
	$i_typ = $file->getColumnIndex("variant_type");
	$i_cod = $file->getColumnIndex("coding_and_splicing");
	
	$parser->log("INSERT/UPDATE VARIANT");
	$c_skip = 0;
	$c_ins = 0;
	$c_upd = 0;
	$t_start = microtime(true);
	
	//insert/update table 'variant'
	$var_ids = array();
	$hash = $db_connect->prepare("INSERT INTO variant (chr, start, end, ref, obs, 1000g, gnomad, gene, variant_type, coding) VALUES (:chr, :start, :end, :ref, :obs, :1000g, :gnomad, :gene, :variant_type, :coding) ON DUPLICATE KEY UPDATE id=id");
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);
		
		//skip variants with too high AF
		if ($max_af>0 && ((is_numeric($row[$i_10g]) && $row[$i_10g]>$max_af) || (is_numeric($row[$i_gno]) && $row[$i_gno]>$max_af)))
		{
			++$c_skip;
			$var_ids[] = -1;
			continue;
		}
		
		//determine variant ID
		$variant_id = $db_connect->getValue("SELECT id FROM variant WHERE chr='".$row[$i_chr]."' AND start='".$row[$i_sta]."' AND end='".$row[$i_end]."' AND ref='".$row[$i_ref]."' AND obs='".$row[$i_obs]."'", -1);
		if ($variant_id!=-1) //update (common case)
		{
			//compose array of meta data
			$metadata = array(
							"1000g" => $row[$i_10g],
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
			$db_connect->bind($hash, "1000g", $row[$i_10g], array(""));
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
	
	if ($max_af>0)
	{
		$parser->log(" skipped variants: $c_skip (AF>$max_af)");
	}
	$parser->log(" inserted variants: $c_ins");
	$parser->log(" updated variants: $c_upd");
	$parser->log(" done ".time_readable(microtime(true) - $t_start));
	
	return $var_ids;
}

//extract variant quality (snp_q) from the GSvar 'quality' column
function getVariantQuality($row, $i_qual)
{
	$quality_parts = explode(";", $row[$i_qual]);
	foreach($quality_parts as $quality_part)
	{
		if (starts_with($quality_part, "QUAL="))
		{
			return explode("=", $quality_part)[1];
		}
	}
	
	trigger_error("Could not extract variant quality from quality entry '".$row[$i_qual]."' of variant ".$row[0].":".$row[1]." ".$row[3].">".$row[4], E_USER_ERROR);
}
		

//check processed sample ID format
$samples = array();
if(preg_match("/^([A-Za-z0-9\-]{4,})_(\d{2})$/", $id, $matches))
{
	if ($mode=="germline")
	{
		$samples["germline"] = array("name"=>($matches[1]."_".$matches[2]), "sname"=>$matches[1],"pname"=>(int)$matches[2]);
	}
	else
	{
		$samples["tumor"] = array("name"=>($matches[1]."_".$matches[2]), "sname"=>$matches[1],"pname"=>(int)$matches[2]);
	}
}
else if(preg_match("/^([A-Za-z0-9]{4,})_(\d{2})-([A-Za-z0-9]{4,})_(\d{2})$/", $id, $matches) && $mode=="somatic")
{
	$samples["tumor"] = array("name"=>($matches[1]."_".$matches[2]),"sname"=>$matches[1],"pname"=>(int)$matches[2]);
	$samples["germline"] = array("name"=>($matches[3]."_".$matches[4]),"sname"=>$matches[3],"pname"=>(int)$matches[4]);
}
else trigger_error("'$id' is not a valid processed sample ID!", E_USER_ERROR);

//load variant file
$file = Matrix::fromTSV($var);

//database connection
$db_connect = DB::getInstance($db);

//convert processed sample ID to NGSD ids
$hash = $db_connect->prepare("SELECT s.id as sid, ps.id as psid FROM sample as s, processed_sample as ps WHERE s.id = ps.sample_id and s.name = :sname and ps.process_id = :pname;");
foreach($samples as $key => $sample)
{
	$db_connect->bind($hash, 'sname', $sample["sname"]);
	$db_connect->bind($hash, 'pname', $sample["pname"]);
	$db_connect->execute($hash, true); 
	$result = $db_connect->fetch($hash);
	if(count($result)!=1)
	{
		trigger_error("Processed sample '".$samples[$key]["name"]."' not found in NGSD!", E_USER_ERROR);
	}
	$samples[$key]["sid"] = $result[0]['sid'];
	$samples[$key]["pid"] = $result[0]['psid'];
}
$db_connect->unsetStmt($hash);

//germline
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

	
	//remove old variants (and store class4/5 variants for check)
	$parser->log("DELETING DETECTED VARIANTS");
	$t_start = microtime(true);
	
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
	$parser->log(" done ".time_readable(microtime(true) - $t_start));

	//abort if there are no variants in the input file
	if ($file->rows()==0)
	{
		$parser->log("No variants imported (empty variant TSV file).");
		exit();
	}
	
	//add missing variants
	$var_ids = updateVariantTable($parser, $db_connect, $file, $max_af);
	
	//insert variants into table 'detected_variant'
	$parser->log("INSERT DETECTED VARIANT");
	$t_start = microtime(true);
	
	$db_connect->beginTransaction();
	$hash = $db_connect->prepare("INSERT INTO detected_variant (processed_sample_id, variant_id, genotype) VALUES ('$psid', :variant_id, :genotype);");
	$i_geno = 5; //genotype is always the 5th column
	$i_typ = $file->getColumnIndex("variant_type");
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);
		$variant_id = $var_ids[$i];
		
		//skip high-AF variants
		if ($variant_id==-1) continue;

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

//somatic - tumor-normal pair
else if($mode=="somatic" && isset($samples["germline"]))
{
	// check if variants were already imported for this PID
	$count_old = $db_connect->getValue("SELECT count(id) FROM detected_somatic_variant WHERE processed_sample_id_tumor='".$samples["tumor"]["pid"]."' AND processed_sample_id_normal='".$samples["germline"]["pid"]."'");
	
	//throw error if not forced to overwrite and terms already exist
	if($count_old!=0 && !$force)
	{
		trigger_error("Variants were already imported for PIDs '".$samples["tumor"]["name"]."-".$samples["germline"]["name"]."'. Use the flag '-force' to overwrite them.", E_USER_ERROR);
	}
	$parser->log("Found $count_old variants for PID '".$samples["tumor"]["name"]."-".$samples["germline"]["name"]."' already in DB.");

	//remove old variants
	if ($count_old!=0 && $force)
	{
		//remove old variants
		$db_connect->executeStmt("DELETE FROM detected_somatic_variant WHERE processed_sample_id_tumor='".$samples["tumor"]["pid"]."' AND processed_sample_id_normal='".$samples["germline"]["pid"]."'");
		$parser->log("Deleted old variants because the flag '-force' was used.");
	}

	//abort if there are no variants in the input file
	if ($file->rows()==0)
	{
		$parser->log("No variants imported (empty variant TSV file).");
		exit();
	}

	//add missing variants
	$var_ids = updateVariantTable($parser, $db_connect, $file);
	
	//insert variants into table 'detected_somatic_variant'
	$i_frq = $file->getColumnIndex("tumor_af");
	$i_dep = $file->getColumnIndex("tumor_dp");
	$i_qual = $file->getColumnIndex("quality");
	$hash = $db_connect->prepare("INSERT INTO detected_somatic_variant (processed_sample_id_tumor,processed_sample_id_normal,variant_id,variant_frequency,depth,quality_snp) VALUES ('".$samples["tumor"]["pid"]."','".$samples["germline"]["pid"]."',:variant_id,:variant_frequency,:depth,:snp_q);");
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);
		
		$db_connect->bind($hash, "variant_id", $var_ids[$i]);
		$db_connect->bind($hash, "variant_frequency", $row[$i_frq], array("n/a"));
		$db_connect->bind($hash, "depth", $row[$i_dep]);
		$db_connect->bind($hash, "snp_q", getVariantQuality($row, $i_qual), array("."));
		
		$db_connect->execute($hash, true);
	}
}

//somatic - only tumor
else if($mode=="somatic" && !isset($samples["germline"]))
{
	// check if variants were already imported for this PID
	$count_old = $db_connect->getValue("SELECT count(id) FROM detected_somatic_variant WHERE processed_sample_id_tumor='".$samples["tumor"]["pid"]."' AND processed_sample_id_normal IS NULL");
	
	//throw error if not forced to overwrite and terms already exist
	if($count_old!=0 && !$force)
	{
		trigger_error("Variants were already imported for processed sample '".$samples["tumor"]["name"]."'. Use the flag '-force' to overwrite them.", E_USER_ERROR);
	}
	$parser->log("Found $count_old variants for processed sample '".$samples["tumor"]["name"]."' already in DB.");

	//remove old variants
	if ($count_old!=0 && $force)
	{
		$db_connect->executeStmt("DELETE FROM detected_somatic_variant WHERE processed_sample_id_tumor='".$samples["tumor"]["pid"]."' AND processed_sample_id_normal IS NULL");
		$parser->log("Deleted old variants because the flag '-force' was used.");
	}
	
	//abort if there are no variants in the input file
	if ($file->rows()==0)
	{
		$parser->log("No variants imported (empty variant TSV file).");
		exit();
	}
	
	//add missing variants
	$var_ids = updateVariantTable($parser, $db_connect, $file);

	//insert variants into table 'detected_somatic_variant'
	$i_frq = $file->getColumnIndex("tumor_af");
	$i_dep = $file->getColumnIndex("tumor_dp");
	$i_qual = $file->getColumnIndex("quality");
	$hash = $db_connect->prepare("INSERT INTO detected_somatic_variant (processed_sample_id_tumor,processed_sample_id_normal,variant_id,variant_frequency,depth,quality_snp) VALUES ('".$samples["tumor"]["pid"]."',NULL,:variant_id,:variant_frequency,:depth,:snp_q);");
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);
		
		$db_connect->bind($hash, "variant_id", $var_ids[$i]);
		$db_connect->bind($hash, "variant_frequency", $row[$i_frq], array("n/a"));
		$db_connect->bind($hash, "depth", $row[$i_dep]);
		$db_connect->bind($hash, "snp_q", getVariantQuality($row, $i_qual), array("."));
		
		$db_connect->execute($hash, true);
	}
}



?>
