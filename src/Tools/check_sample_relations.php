<?php
/** 
	@page db_import_sample_relations
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("check_sample_relations", "Checks same sample relationships by verifying genotypes.");
$parser->addString("project", "Samples in this project are checked.", false);
//optional
$parser->addString("roi", "Target region to use.", true, "/mnt/share/data/enrichment/ssSC_v5_2019_07_09.bed");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

$min_correlation = 0.9;

$db = DB::getInstance($db);

$query = <<<SQL
SELECT
CONCAT(s1.name, "_", LPAD(ps1.process_id, 2, "0")) as psamplename1,
CONCAT(s2.name, "_", LPAD(ps2.process_id, 2, "0")) as psamplename2
FROM sample_relations sr
LEFT JOIN sample s1 on s1.id=sample1_id
LEFT JOIN sample s2 on s2.id=sample2_id
LEFT JOIN processed_sample ps1 on ps1.sample_id=s1.id
LEFT JOIN processed_sample ps2 on ps2.sample_id=s2.id
LEFT JOIN project pr1 on pr1.id=ps1.project_id
LEFT JOIN project pr2 on pr2.id=ps2.project_id
LEFT JOIN merged_processed_samples mps1 on mps1.processed_sample_id=ps1.id
LEFT JOIN merged_processed_samples mps2 on mps2.processed_sample_id=ps2.id
WHERE mps1.merged_into IS NULL AND mps2.merged_into IS NULL AND (pr1.name=:project OR pr2.name=:project)
SQL;

$res = $db->executeQuery($query, array('project'=>$project));
$pairs = [];
$checked_samples_in_project = [];
foreach ($res as $row)
{
    $key1 = $row['psamplename1'] . "-" . $row['psamplename2'];
    $key2 = $row['psamplename2'] . "-" . $row['psamplename1'];
    if (in_array($key1, $pairs) || in_array($key2, $pairs))
    {
        print("Skipping pair {$key1}...\n");
        continue;
    }
    array_push($pairs, $key1);
    // run sample check
    // print("Run sample check for pair {$key1}...\n");
    $info1 = get_processed_sample_info($db, $row['psamplename1']);
    $info2 = get_processed_sample_info($db, $row['psamplename2']);
    $bam1 = $info1['ps_bam'];
    $bam2 = $info2['ps_bam'];
    if (!file_exists($bam1))
    {
        trigger_error("Missing BAM file '{$bam1}'!", E_USER_WARNING);
        continue;
    }
    if (!file_exists($bam2))
    {
        trigger_error("Missing BAM file '{$bam2}'!", E_USER_WARNING);
        continue;
    }

    if ($info1['project_name'] === $project)
    {
        array_push($checked_samples_in_project, $row['psamplename1']);
    }
    else if ($info2['project_name'] === $project)
    {
        array_push($checked_samples_in_project, $row['psamplename2']);
    }

    $args_similarity = [
        "-in", $bam1, $bam2,
        "-mode", "bam",
        "-min_cov", 10
    ];
    if (!empty($roi))
    {
        $args_similarity[] = "-roi";
        $args_similarity[] = $roi;
    }
    $output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", implode(" ", $args_similarity), true);
    $n_snps = explode("\t", $output[0][1])[2];
    $correlation = explode("\t", $output[0][1])[3];
    $flag = $correlation < $min_correlation ? "ERROR" : "PASS";
    print(implode("\t", [$row['psamplename1'], $row['psamplename2'], $correlation, $n_snps, $flag]) . "\n");
}

$query = <<<SQL
SELECT
CONCAT(s.name, "_", LPAD(ps.process_id, 2, "0")) as psamplename
FROM processed_sample ps
LEFT JOIN sample s on s.id=ps.sample_id
LEFT JOIN project pr on pr.id=ps.project_id
LEFT JOIN merged_processed_samples mps on mps.processed_sample_id=ps.id
WHERE mps.merged_into IS NULL AND pr.name=:project
SQL;
$res = $db->executeQuery($query, array('project'=>$project));
$samples_in_project = array_column($res, "psamplename");
$checked_samples_in_project = array_unique($checked_samples_in_project);
$unchecked_samples_in_project = array_diff($samples_in_project, $checked_samples_in_project);
print("\nSamples checked:\n" . implode("\n", $checked_samples_in_project) . "\n");
print("\nSamples in project but no relation found:\n" . implode("\n", $unchecked_samples_in_project) . "\n");

?>