<?php
/** 
	@page db_run_status
*/
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("db_run_status", "\$Rev: 592 $", "Prints status information about run folders in the /mnt/raw-data/");
extract($parser->parse($argv));

function find_run($fcid)
{	
	$db_connect = DB::getInstance("NGSD");
	$hash = $db_connect->prepare("SELECT name, quality, status, comment FROM sequencing_run WHERE fcid=:fcid");
	$db_connect->bind($hash, "fcid", $fcid);
	$db_connect->execute($hash, TRUE);
	return $db_connect->fetch($hash);
}


//find/validate backups
$backups = array();
$tmp = array();
exec("find /mnt/archive/runs* -maxdepth 1 -name \"*.tgz\" -or -name \"*.tar.gz\"", $tmp);
foreach($tmp as $t)
{
	$mb = filesize($t) / 1024^2;
	if($mb<500)
	{
	
		print "#Warning: $t is too small ($mb MB)!\n";
		continue;
	}
	$backups[] = basename($t);
}

//init
print "#folder\tfcid\tname\tquality\tstatus\tbackup\tcomment\n";
$folders = array();
exec("find /mnt/raw-data/ -maxdepth 2 -type d -name \"[1]*\"", $folders);
foreach($folders as $folder)
{
	//tokenize run folder name
	$parts = explode("_", basename($folder));
	$c = count($parts);

	//check if run number was appended
	$fcid = $parts[$c-1];
	if (preg_match("/^[0-9][0-9][0-9][0-9][0-9]$/", $fcid))
	{
		$fcid = $parts[$c-2];
	}

	//remove prepended 'A'/'B' from HiSeq
	if (strlen($fcid)>5 && ($fcid[0]=='A' || $fcid[0]=='B'))
	{
		$fcid = substr($fcid, 1);
	}
	
	//remove '000000000-'
	$fcid = str_replace("000000000-", "", $fcid);

	//determine run
	$run = find_run($fcid);
	if ($run!=array())
	{
		$name = str_replace("#", "", $run[0]["name"]);
		$quality = $run[0]["quality"];
		$status = $run[0]["status"];
		$comment = strtr($run[0]["comment"], array("\n"=>"", "\r"=>""));
	}
	else
	{
		$name = "";
		$quality = "";
		$status = "";
		$comment = "";
	}

	//check if backup was performed
	$backup="backup_missing";
	foreach($backups as $b)
	{
		if (($name!="" && ends_with($b, $name.".tar.gz")) || ends_with($b, $fcid.".tar.gz"))
		{
			$backup = "backup_done";
			break;
		}
	}

	print "$folder\t$fcid\t$name\t$quality\t$status\t$backup\t$comment\n";
}

?>