<?php
/** 
	@page export_samplesheet
	
	@todo 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_samplesheet", "Exports a bcl2fastq2-compatible sample sheet.");
$parser->addString("run", "Run name.", false);
$parser->addOutfile("out", "Output file in CSV format.", false);
$parser->addString("lanes", "Comma-separated list of lane numbers to use.", true, "1,2,3,4,5,6,7,8");
$parser->addInt("mid1_len", "Numer of bases to use from MID 1.", true, -1);
$parser->addInt("mid2_len", "Numer of bases to use from MID 2.", true, -1);
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db = DB::getInstance($db);
$lanes = explode(",", $lanes);
$output = array();

//execute query
$res = $db->executeQuery("SELECT CONCAT(s.name, '_', LPAD(ps.process_id,2,'0')) as psname, (select sequence from mid where id = ps.mid1_i7) as mid1_i7, (select sequence from mid where id = ps.mid2_i5) as mid2_i5, ps.comment as pscomment,ps.lane as lane, s.comment as scomment, p.name as pname FROM processed_sample as ps, sample as s, genome as g, project as p, sequencing_run as sr, device as d, processing_system as sys WHERE ps.sample_id = s.id and ps.project_id = p.id and sys.id = ps.processing_system_id and sys.genome_id = g.id and sr.id = ps.sequencing_run_id AND sr.device_id = d.id and sr.name = :run ORDER BY psname", array("run" => $run));

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
		$value = str_ireplace(array("\"","ö","ä","ß","ü"),array("","oe","ae","ss","ue"), $value);
		//trim
		$value = trim($value);
		
		$res[$index][$key] = $value;
	}
	
	//convert index 2 to reverse-complement
	$res[$index]["mid2_i5"] = rev_comp($res[$index]["mid2_i5"]);
}

//generate sample sheet
$barcodes = array();
$output[] = "[Data]";
$output[] = "Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description";
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
		$barcodes[$lane][] = array($mid1, $mid2, $name);
		$output[] = $lane.",Sample_".$row["psname"].",".$name.",,,,".$mid1.",,".$mid2.",".$row["pname"].",".$row["pscomment"]." ".$row["scomment"];
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
			if ($dist==0)
			{
				trigger_error("Barcode clash on lane {$lane}: {$name_i} ({$mid1_i},{$mid2_i}) / {$name_j} ({$mid1_j},{$mid2_j})", E_USER_ERROR);
			}
			else if ($dist<=2)
			{
				trigger_error("Low barcode distiance ({$dist}) on lane {$lane}: {$name_i} ({$mid1_i},{$mid2_i}) / {$name_j} ({$mid1_j},{$mid2_j})", E_USER_WARNING);
			}
		}
	}
}

?>