<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("transfer_run_mibi", "Transfers raw run data to Mibi.");
$parser->addString("run", "Run name.", false);
$parser->addInfile("in",  "Run folder.", false);
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addFlag("dry_run", "Execute rsync in dry-run mode.");

extract($parser->parse($argv));

//trim trailing slash, important for rsync
$in = rtrim($in, "/");

//write IMGAG_RunInfo.txt file
$db = DB::getInstance($db);
$res = $db->executeQuery("SELECT name, fcid, start_date, end_date, pool_molarity, comment FROM sequencing_run WHERE name=:name", ['name' => $run]);
$lines = [];
foreach ($res[0] as $k => $v)
{
    $lines[] = "{$k}: {$v}";
}
file_put_contents("{$in}/IMGAG_RunInfo.txt", implode("\n", $lines) . "\n");

$target = rtrim(get_path("mibi_target"), "/");

//copy key to ensure restrictive permissions required for SSH
$keyfile = get_path("mibi_key");
$key_temp = $parser->tempFile();
touch($key_temp);
chmod($key_temp, 0600);
file_put_contents($key_temp, file_get_contents($keyfile));

$args = [
    "--recursive",
    "--verbose",
    "--human-readable",
    "--rsh 'ssh -i {$key_temp}'",
    $dry_run ? "--dry-run" : "",
    $in,
    $target
];
$parser->exec("rsync", implode(" ", $args));

?>