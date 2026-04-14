#! /usr/bin/php
<?php 

/** 
	@page backup_restore
*/


require_once("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("backup_restore", "Restores a file from TSM .");
$parser->addString("file",  "Backup log file or name of tar.gz file to restore.", false);
$parser->addEnum("mode",  "Restore mode which type of backup will be restored.", false, ["runs", "users", "projects"]);
$parser->addString("out",  "Folder where the file will be restored to.", true, "/mnt/storage1/raw_data/_restored");

extract($parser->parse($argv));


function query($file)
{
	echo "\n*** The archive contains the following versions (only the most recent one will be restored):\n";
	exec("dsmc query archive {$file}", $output, $res);
	
	if ($res == 0)
	{
		$offset =0;
		for ($i=0; $i<count($output); $i++)
		{
			if (contains($output[$i], "-----"))
			{
				$offset = $i-1;
				break;
			}
		}
		
		echo implode("\n", array_slice($output, $offset))."\n\n";
	}
	else
	{
		trigger_error("Error retrieving file versions. Error code $res.", E_USER_ERROR);
	}
}

function restore($source, $to)
{
	passthru("dsmc retrieve \"{$source}\"  \"{$to}\"", $res);
	
	if ($res == 0)
	{
		echo "Setting permissions for file...\n";
		exec("chgrp \"f_ad_bi_l_medgen_access_storages\" {$to}");
		exec("chmod 664 {$to}");
		echo "File was successfully restored!\n";
		echo "{$to}\n";
	}
	else
	{
		trigger_error("Error restoring the file. Error code: $res.", E_USER_ERROR);
	}
}


# *** MAIN ***
putenv("DSM_DIR=/opt/tsm-archive/config/");
putenv("DSM_CONFIG=/opt/tsm-archive/config/dsm.opt");

#check script permissions are 755:
$scriptPath = get_included_files()[0];
$permissions_readable = substr(sprintf('%o', fileperms($scriptPath)), -4);
if ($permissions_readable != "0755") trigger_error("This script should have the permissions set to 755! Please have an admin correct the permissions before using it.", E_USER_ERROR);
	
#check if user is root
$executingUser = posix_getpwuid(posix_geteuid())["name"];
if ($executingUser != "root")
{
	trigger_error("This script can only be run as root (sudo).", E_USER_ERROR);
}

if(!is_dir($out)) trigger_error("Given out path '$out' does not exist or is not a directory.", E_USER_ERROR);

$filename = basename($file, ".log");
if (! ends_with(".tar.gz", $file)) $filename .= ".tar.gz";

//check filename for forbidden symbols:
$matched = preg_match_all('/[^A-Za-z0-9\_\.\-]/', $filename);
if ($matched != 0)
{
	trigger_error("Filename '{$filename}' contains forbidden symbols only 'A-Z', 'a-z', '0-9', '_', '.' and '-' are allowed.", E_USER_ERROR);
}

$file = "/mnt/SRV018/raw_data_archive/{$mode}/{$filename}";
$to = "{$out}/{$filename}";

# checks:

if (! is_dir($out)) trigger_error("The out folder '$out' does not exist!", E_USER_ERROR);
if (file_exists($to)) trigger_error("The file $filename already exists in the outFolder!", E_USER_ERROR);
if (file_exists("/mnt/storage1/raw_data_archive/runs/{$filename}")) trigger_error("The file still exists in the archiving folder: /mnt/storage1/raw_data_archive/runs/{$filename}", E_USER_ERROR);
$permissions_readable = substr(sprintf('%o', fileperms($out)), -4);
if ($permissions_readable != "0777") trigger_error("The out folder needs permissions set to 777 during the restore. Pleace set the accordingly and restart the script.", E_USER_ERROR);



echo "Searching for $filename in archive.\n";

query($file);

$input = readline("Press ENTER, to start the restoration of the file. Press Strg+C to abort.");

restore($file, $to);

?>
