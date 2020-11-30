<?php
/** 
	@page export_samplesheet
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("generate_demux_details", "Exports Stats.json information for each project and create sample sheet for MultiQC report");
$parser->addString("run", "Run name.", false);
$parser->addFlag("execute_multiqc", "Executes multiqc for each project in the run.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db = DB::getInstance($db);
$output = array();

//execute query
$res = $db->executeQuery("SELECT CONCAT(s.name, '_', LPAD(ps.process_id,2,'0')) as psname, (select sequence from mid where id = ps.mid1_i7) as mid1_i7, (select sequence from mid where id = ps.mid2_i5) as mid2_i5, ps.comment as pscomment,ps.lane as lane, s.comment as scomment, p.name as pname, p.type as ptype, s.name_external as ename, s.name as sname, s.sample_type as stype, s.species_id as sspecies FROM processed_sample as ps, sample as s, genome as g, project as p, sequencing_run as sr, device as d, processing_system as sys WHERE ps.sample_id = s.id and ps.project_id = p.id and sys.id = ps.processing_system_id and sys.genome_id = g.id and sr.id = ps.sequencing_run_id AND sr.device_id = d.id and sr.name = :run ORDER BY psname", array("run" => $run));

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
}

//generate sample sheets
$projects = array();
$project_locations = array();
foreach($res as $row)
{
	if(array_key_exists($row['pname'], $projects)) {
		$projects[$row['pname']][] = $row['psname']."\t".$row['ename']."\t".$row['stype'];
	} else {
		$projects[$row['pname']][] = "Internal Name\tExternal Name\tType";
		$projects[$row['pname']][] = $row['psname']."\t".$row['ename']."\t".$row['stype'];
		$project_locations[$row['pname']] = "/mnt/projects/".$row['ptype']."/".$row['pname'];
	}
}

//store sample sheets
foreach($projects as $key =>$value)
{
	file_put_contents($project_locations[$key]."/".$key."_samplesheet.tsv", implode("\n", $value));
}

//generate Stats.json files
$projects = array();
foreach($res as $row)
{
	if(array_key_exists($row['pname'], $projects)) {
		$projects[$row['pname']][] = $row['psname'];
	} else {
		$projects[$row['pname']][] = array();
		$projects[$row['pname']][] = $row['psname'];
	}
}

$demux_file = 'Unaligned/Stats/Stats.json';
$unset_keys = array();

if(file_exists($demux_file)) {
    $file = file_get_contents('Unaligned/Stats/Stats.json');
	foreach($projects as $project => $samples)
	{
		$data = json_decode($file);
		foreach($data->ConversionResults as $conversion_result) {
			foreach($conversion_result->DemuxResults as $key => $value) {
				foreach($value as $k => $v) {
					if($k == "SampleName") {
						if(!in_array($v,$samples)) {
							$unset_keys[]=$key;
						}
					}
				}
				$value->SampleId = $value->SampleName;
			}
			
			foreach($unset_keys as $u_key) {
				unset($conversion_result->DemuxResults[$u_key]);
			}
			$conversion_result->DemuxResults = array_values($conversion_result->DemuxResults);
		}
		
		//store Stats.json files
		file_put_contents($project_locations[$project]."/Stats.json", json_encode($data));
	}
}

if ($execute_multiqc)
{
	foreach($projects as $key =>$value)
	{
		exec("multiqc ".$project_locations[$key]."/ --sample-names ".$project_locations[$key]."/".$key."_samplesheet.tsv -f -o ".$project_locations[$key]);
	}
}

?>
