<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("2021_04_subpanel_import", "Imports subpanels.");
$parser->addString("folder", "Subpanel base folder.", false);
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db = DB::getInstance($db);

$user2id = [];
$res = $db->executeQuery("SELECT id, user_id FROM user");
foreach($res as $row)
{
	$user2id[strtolower($row['user_id'])] = $row['id'];
}

//import
function import($folder, $is_archive)
{
	global $db;
	global $user2id;
		
	list($files) = exec2("ls $folder");
	foreach($files as $file)
	{
		if (!ends_with($file, ".bed")) continue;
		if (ends_with($file, "_amplicons.bed")) continue;
		
		$name = trim(substr($file, 0, -4));
		$roi_file = $folder."/".$file;
		$gene_file = $folder."/".substr($file, 0, -4)."_genes.txt";
		if (!file_exists($gene_file))
		{
			print "$file: genes file missing!\n";
			continue;
		}
		
		//parse name
		while(contains($file, "__")) $file = strtr($file, ["__"=>"_"]);
		$parts = explode("_", substr($file, 0, -4));
		if (count($parts)<4)
		{
			"$file: too few parts\n";
			continue;
		}
		list($mode_ext, $user, $date) = array_slice($parts, -3, 3); 
		
		$user = strtolower($user);
		if (!isset($user2id[$user]))
		{
			print "$file: unknown user $user!\n";
			continue;
		}
		$user_id = $user2id[$user];

		$matches = [];
		if (!preg_match("/([a-z]+)([0-9]+)/", $mode_ext, $matches))
		{
			print "$file: could not split mode/extension!\n";
			continue;
		}
		$mode = $matches[1];
		$ext = $matches[2];
		
		$genes = file_get_contents($gene_file);
		$roi = file_get_contents($roi_file);
		
		$db->executeStmt("INSERT INTO `subpanels`(`name`, `created_by`, `created_date`, `mode`, `extend`, `genes`, `roi`, `archived`) VALUES ('$name', $user_id, '$date', '$mode', $ext, '$genes', '$roi', $is_archive) ON DUPLICATE KEY UPDATE archived='$is_archive'");
	}
}

print "\n";
print "Importing active subpanels...\n";
import($folder, "0");

print "\n";
print "Importing archived subpanels...\n";
import($folder."/archive/", "1");

?>