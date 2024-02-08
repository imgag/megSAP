<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("rc_annotate_expr", "Annotates a read count table with relative expression values.");
$parser->addString("name", "Processed sample name used to find related samples.", false);
$parser->addInfile("in", "Input gene counts table.", false);
$parser->addOutfile("out", "Output gene counts table.", false);
$parser->addString("cohort", "Specify to save expression values of cohort.", true, "");
$parser->addString("stats", "Specify to save statistics of expression values of cohort.", true, "");
$parser->addString("corr", "Specify to save sample--cohort correlation.", true, "");
$parser->addInfile("in_files", "Tab-separated file with sample name and corresponding path to count file per line.", true);

$parser->addFlag("somatic", "Enable somatic mode: find related samples by ICD10 code or HPO term id.");
$parser->addString("hpa_tissue", "HPA reference tissue.", true, "");
$parser->addString("hpa_corr", "Specify to save sample--reference tissue correlation.", true, "");
$parser->addString("hpa_ref", "HPA reference expression file.", true, get_path("data_folder") . "/dbs/gene_expression/rna_tissue_consensus_v23.tsv");
$parser->addEnum("db", "Database to connect to.", true, db_names(), "NGSD");

extract($parser->parse($argv));


function get_related_samesample(&$db_conn, $sample_id)
{
	$related_samples = [];
	$res = $db_conn->executeQuery("SELECT * FROM sample_relations WHERE relation='same sample' AND (sample1_id=:sid OR sample2_id=:sid)", array("sid" => $sample_id));
	foreach($res as $row)
	{
		$related_samples[] = $row['sample1_id'] != $sample_id ? $row['sample1_id'] : $row['sample2_id'];
	}

	return $related_samples;
}

if (!isset($in_files))
{
	$db_conn = DB::getInstance($db);

	//sample and processed sample values to match
	list($sample_name, $process_id) = explode("_", $name."_");
	$psid = get_processed_sample_id($db_conn, $name);
	$sample_id = $db_conn->getValue("SELECT id FROM sample WHERE name='{$sample_name}'");
	$processing_system_id = $db_conn->getValue("SELECT processing_system_id FROM processed_sample WHERE id={$psid}");
	$project_id = $db_conn->getValue("SELECT project_id FROM processed_sample WHERE id={$psid}");

	//collected samples and count files
	$counts_paths = [];

	if ($somatic)
	{
		//somatic use case:
		//extract HPO term id and ICD10 code of RNA sample, or of samples with "same sample" relation
		$related_samples = get_related_samesample($db_conn, $sample_id);
		$related_samples[] = $sample_id;

		$icd10 = [];
		$hpo = [];
		$query = "SELECT disease_info FROM sample_disease_info WHERE sample_id=:sample_id AND type=:type";
		foreach ($related_samples as $sample_id)
		{
			$res = $db_conn->executeQuery($query, ["sample_id"=>$sample_id, "type"=>"ICD10 code"]);
			$icd10 = array_unique(array_merge($icd10, array_column($res, "disease_info")));

			$res = $db_conn->executeQuery($query, ["sample_id"=>$sample_id, "type"=>"HPO term id"]);
			$hpo = array_unique(array_merge($hpo, array_column($res, "disease_info")));
		}

		if (count($icd10) >= 2) trigger_error("More than one ICD10 code found!", E_USER_WARNING);
		if (count($hpo) >= 2) trigger_error("More than one HPO term found!", E_USER_WARNING);

		if (count($icd10) == 0) trigger_error("No ICD10 code found!", E_USER_ERROR);
		if (count($hpo) == 0) trigger_error("No HPO term found!", E_USER_ERROR);

		//find related samples, matching the following values:
		// trigger_error("processing system id: {$processing_system_id}", E_USER_NOTICE);
		// trigger_error("project id: {$project_id}", E_USER_NOTICE);
		// trigger_error("ICD10 code: {$icd10[0]}", E_USER_NOTICE);
		// trigger_error("HPO term id: {$hpo[0]}", E_USER_NOTICE);

		$sql = <<<SQL
SELECT DISTINCT ps.id, CONCAT(s.name, "_", LPAD(ps.process_id, 2, "0")) as psample
FROM processed_sample ps
LEFT JOIN sample s on ps.sample_id=s.id
LEFT JOIN sample_relations sr ON s.id=sr.sample1_id OR s.id=sr.sample2_id
LEFT JOIN sample_disease_info sdi ON s.id=sdi.sample_id OR sr.sample1_id=sdi.sample_id OR sr.sample2_id=sdi.sample_id
WHERE
	ps.processing_system_id=:psys_id
	AND
	ps.project_id=:project_id
	AND
	(sr.relation="same sample" OR sr.relation IS NULL)
	AND
	((sdi.type="ICD10 code" AND sdi.disease_info=:icd10)
	OR
	(sdi.type="HPO term id" AND sdi.disease_info=:hpo));
SQL;

		$res = $db_conn->executeQuery($sql, ["psys_id"=>$processing_system_id, "project_id"=>$project_id, "icd10"=>$icd10[0], "hpo"=>$hpo[0]]);
		$psamples = array_unique(array_column($res, "psample"));
		foreach ($psamples as $psample)
		{
			$info = get_processed_sample_info($db_conn, $psample);

			//skip bad quality samples
			if ($info['ps_quality']=="bad") continue;

			$counts_path = "{$info['ps_folder']}{$psample}_counts.tsv";
			if (file_exists($counts_path))
			{
				$counts_paths[$psample] = $counts_path;
			}
			else
			{
				trigger_error("Missing count file: '{$counts_path}'", E_USER_WARNING);
			}
		}
	}
	else
	{
		//default use case:
		//find related samples, matching the following values:
		trigger_error("processing system id: {$processing_system_id}", E_USER_NOTICE);
		trigger_error("project id: {$project_id}", E_USER_NOTICE);

		$sql = <<<SQL
SELECT DISTINCT ps.id, CONCAT(s.name, "_", LPAD(ps.process_id, 2, "0")) as psample
FROM processed_sample ps
LEFT JOIN sample s on ps.sample_id=s.id
WHERE
	ps.processing_system_id=:psys_id
	AND
	ps.project_id=:project_id
	AND
	ps.quality!='bad'
SQL;

		$res = $db_conn->executeQuery($sql, ["psys_id"=>$processing_system_id, "project_id"=>$project_id]);
		$psamples = array_unique(array_column($res, "psample"));
		foreach ($psamples as $psample)
		{
			$info = get_processed_sample_info($db_conn, $psample);
			$counts_path = "{$info['ps_folder']}{$psample}_counts.tsv";
			if (file_exists($counts_path))
			{
				$counts_paths[$psample] = $counts_path;
			}
			else
			{
				trigger_error("Missing count file: '{$counts_path}'", E_USER_WARNING);
			}
		}
	}


	//write file of file names
	$lines = [];
	ksort($counts_paths);
	foreach ($counts_paths as $psample => $path)
	{
		$lines[] = "{$psample}\t{$path}\n";
	}
	$fofn_tmp = $parser->tempFile("counts.fofn");
	file_put_contents($fofn_tmp, $lines);
	$in_files = $fofn_tmp;
}

//annotate cohort information
$args = [
	"cohort",
	"--counts", $in,
	"--counts_out", $out,
	"--samples", $in_files
];
if ($cohort !== "") $args[] = "--cohort {$cohort}";
if ($stats !== "") $args[] = "--stats {$stats}";
if ($corr !== "") $args[] = "--corr {$corr}";
$parser->exec(get_path("python3")." ".repository_basedir()."/src/NGS/rc_calc_expr.py", implode(" ", $args), true);

//TODO: use ngs-bits tool

//annotate HPA information
if ($hpa_tissue !== "")
{
	$args = [
		"hpa",
		"--counts", $out,
		"--counts_out", $out,
		"--hpa", $hpa_ref,
		"--prefix", "hpa_",
		"--tissue", "'{$hpa_tissue}'"
	];
	if ($hpa_corr !== "") $args[] = "--corr {$hpa_corr}";
	$parser->exec(get_path("python3")." ".repository_basedir()."/src/NGS/rc_calc_expr.py", implode(" ", $args), true);
}
