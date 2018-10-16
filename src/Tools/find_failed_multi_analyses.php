<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

list($stdout) = exec2("find /mnt/projects/ -maxdepth 3 -type d -name \"Trio_*\" -or -name \"Multi_*\"");
foreach($stdout as $folder)
{
	$gsvar_file = "{$folder}/".strtolower(explode("_", basename($folder))[0]).".GSvar";
	if (!file_exists($gsvar_file))
	{
		$files = glob("$folder/*.*");
		print "{$folder}\n";
		foreach($files as $file)
		{
			print "    {$file}\n";
		}
	}
}


?>
