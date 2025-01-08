<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("rename_sample", "Renames sample (file system and database).");
$parser->addString("old", "Old sample name.", false);
$parser->addString("new", "New sample name.", false);
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

$db = DB::getInstance($db);

//check sample name
$res = $db->executeQuery("SELECT s.* FROM sample as s WHERE name=:name", array("name" => $old));
if(count($res) === 0)
{
	trigger_error("No sample found for name {$old}!", E_USER_ERROR);
}
if(count($res) >= 2)
{
	trigger_error("Multiple samples found for name {$old}!", E_USER_ERROR);
}
$res_new = $db->executeQuery("SELECT s.* FROM sample as s WHERE name=:name", array("name" => $new));
if(count($res_new) !== 0)
{
	trigger_error("Sample found for new name {$new}!", E_USER_ERROR);
}

//move files
$processed_samples = $db->executeQuery("SELECT ps.id, ps.process_id, CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as psample FROM processed_sample as ps, sample as s WHERE ps.sample_id=:id AND s.id=:id AND ps.id NOT IN (SELECT processed_sample_id FROM merged_processed_samples)", array("id" => $res[0]["id"]));
foreach ($processed_samples as $ps)
{
	$old_ps = $ps["psample"];
	$info_old = get_processed_sample_info($db, $old_ps);
	$dir = $info_old["ps_folder"];

	$fastq_files = glob("${dir}/{$old}*.fastq.gz");
	foreach ($fastq_files as $file)
	{
		$basename_new = preg_replace("/^{$old}/", $new, basename($file));
		$path_new = dirname($file) . "/" . $basename_new;
		$parser->exec("mv", "{$file} {$path_new}");
	}

	$old_files = glob("${dir}/{$old}*");
	$old_dir = "${dir}/+rename";
	create_directory($old_dir);
	foreach ($old_files as $file)
	{
		$parser->exec("mv", "{$file} ${old_dir}/");
	}

	$dir_new = dirname($dir) . "/" . preg_replace("/^Sample_{$old}/", "Sample_{$new}", basename($dir));
	$parser->exec("mv", "{$dir} {$dir_new}");
}

//update DB, append old name in external name
$new_external = implode(", ", [ $res[0]["name_external"], $old ]);
$db->executeStmt("UPDATE sample SET name=:new, name_external=:new_external WHERE id=:id", array("new" => $new, "new_external" => $new_external, "id" => $res[0]["id"]));

?>