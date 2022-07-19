<?php
/** 
	@page export_samplesheet
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_samplesheet", "Exports a bcl2fastq2-compatible sample sheet.");
$parser->addString("run", "Run name.", false);
$parser->addOutfile("out", "Output file in CSV format.", false);
$parser->addString("lanes", "Comma-separated list of lane numbers to use.", true, "1,2,3,4,5,6,7,8");
$parser->addInt("mid1_len", "Number of bases to use from MID 1.", true, -1);
$parser->addInt("mid2_len", "Number of bases to use from MID 2.", true, -1);
$parser->addFlag("mid2_no_rc", "Disables reverse-complement of MID 2.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db = DB::getInstance($db);
$lanes = explode(",", $lanes);
$output = array();

//execute query
$res = $db->executeQuery(<<<SQL
SELECT
	CONCAT(s.name, '_', LPAD(ps.process_id,2,'0')) as psname,
	(select sequence from mid where id = ps.mid1_i7) as mid1_i7,
	(select sequence from mid where id = ps.mid2_i5) as mid2_i5,
	(select name from mid where id = ps.mid1_i7) as mid1_i7_name,
	(select name from mid where id = ps.mid2_i5) as mid2_i5_name,
	ps.comment as pscomment,
	ps.lane as lane,
	s.comment as scomment,
	p.name as pname
FROM
	processed_sample as ps,
	sample as s,
	project as p,
	sequencing_run as sr
WHERE
	ps.sample_id = s.id and
	ps.project_id = p.id and
	sr.id = ps.sequencing_run_id and
	sr.name = :run
ORDER BY psname
SQL
, array("run" => $run));

if(count($res)==0)
{
	trigger_error("No processed samples found for run {$run}!", E_USER_ERROR);
}

//clean data
foreach($res as $index => $row)
{
	
	//sanitize text data
	foreach($row as $key => $value)
	{
		//merge whitespaces
		$value = preg_replace('/\s+/',' ', $value);
		//replace comma by dot
		$value = strtr($value, ",", ".");
		//special characters
		$value = str_ireplace(array("\"","ö","ä","ß","ü","µ"),array("","oe","ae","ss","ue","u"), $value);
		//trim
		$value = trim($value);
		
		$res[$index][$key] = $value;
	}
	
	//convert index 2 to reverse-complement
	if (!$mid2_no_rc)
	{
		$res[$index]["mid2_i5"] = rev_comp($res[$index]["mid2_i5"]);
	}
}

//generate sample sheet
$barcodes = array();
$output[] = "[Data]";
$output[] = "Lane,Sample_ID,Sample_Name,Sample_Project,index,index2,comment,index_name,index2_name";
foreach($res as $row)
{
	$lanes_sample = explode(".", $row['lane']);
	foreach($lanes_sample as $lane)
	{
		if (!in_array($lane, $lanes)) continue;
		
		$mid1 = $row["mid1_i7"];
		if ($mid1_len!=-1)
		{
			$mid1 = substr($mid1, 0, $mid1_len);
		}
		$mid2 = $row["mid2_i5"];
		if ($mid2_len!=-1)
		{
			$mid2 = substr($mid2, 0, $mid2_len);
		}
		$name = $row["psname"];

		$mids_i7 = [ $mid1 ];
		$mids_i7_names = [ $row["mid1_i7_name"] ];
		if (starts_with($row["pscomment"], "add_mids:"))
		{
			$add_mids = explode(".", substr($row["pscomment"], strpos($row["pscomment"], ":")+1));
			foreach ($add_mids as $mid)
			{
				$mid_trimmed = trim($mid);
				$res = $db->executeQuery("SELECT name, sequence FROM mid WHERE name = :name", array("name" => $mid_trimmed));
				if (count($res) === 0)
				{
					trigger_error("Invalid MID name: '{$mid_trimmed}'", E_USER_WARNING);
				}
				else
				{
					$mids_i7[] = $res[0]["sequence"];
					$mids_i7_names[] = $res[0]["name"];
				}
			}
		}

		foreach (array_combine($mids_i7_names, $mids_i7) as $mid1_name => $mid1)
		{
			$barcodes[$lane][] = array($mid1, $mid2, $name);
			$output[] = implode(',', [
				$lane,
				"Sample_{$name}",
				$name,
				$row["pname"],
				$mid1,
				$mid2,
				trim(implode(" ", [ $row["pscomment"], $row["scomment"] ])),
				$mid1_name,
				$row["mid2_i5_name"]
			]);
		}
	}
}

//store sample sheet
file_put_contents($out, implode("\n", $output));

//barcode clash check
foreach($barcodes as $lane => $data)
{
	for($i=0; $i<count($data); ++$i)
	{
		for($j=$i+1; $j<count($data); ++$j)
		{
			list($mid1_i, $mid2_i, $name_i) = $data[$i];
			list($mid1_j, $mid2_j, $name_j) = $data[$j];
			$dist1 = levenshtein($mid1_i, $mid1_j);
			$dist2 = levenshtein($mid2_i, $mid2_j);
			
			$dist = $dist1 + $dist2;
			
			//Check whether dual and single barcodes are mixed (if so, the distance is only the first barcode for one of the samples)
			if((empty($mid2_i) && !empty($mid2_j)) || (!empty($mid2_i) && empty($mid2_j))) $dist = $dist1;
			
			if ($dist==0)
			{
				trigger_error("Barcode clash on lane {$lane}: {$name_i} ({$mid1_i},{$mid2_i}) / {$name_j} ({$mid1_j},{$mid2_j})", E_USER_ERROR);
			}
			else if ($dist<=2)
			{
				trigger_error("Low barcode distance ({$dist}) on lane {$lane}: {$name_i} ({$mid1_i},{$mid2_i}) / {$name_j} ({$mid1_j},{$mid2_j})", E_USER_WARNING);
			}
		}
	}
}

//check recipe and suggest basemask and make UMI warning
$recipe = $db->getValue("SELECT recipe FROM sequencing_run WHERE name='" . $run . "'");
$cycles = array_map("intval", array_map("trim", explode("+", $recipe)));
if (count($cycles) === 4)
{
	list($cycles_R1, $cycles_I1, $cycles_I2, $cycles_R2) = $cycles;
	if ($cycles_I1 > $cycles_I2)
	{
		$additional_cycles = $cycles_I1 - $cycles_I2;
		print(
			"Recipe: {$recipe}\n" .
			"More cycles in index 1 than index 2! Please check for UMI in index 1 and use correct base mask.\n" .
			"Suggestion for base mask: Y{$cycles_R1},I{$cycles_I2}Y{$additional_cycles},I{$cycles_I2},Y{$cycles_R2}\n"
		);
	}
}

?>
