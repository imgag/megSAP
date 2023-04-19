#! /usr/bin/php
<?php 

/** 
	@page backup_restore
*/


require_once("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("backup_restore", "Restores a file from TSM .");
$parser->addString("file",  "Backup log file or name of file to restore.", false);
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
	exec("dsmc retrieve \"{$source}\"  \"{$to}\"", $output, $res);
	
	if ($res == 0)
	{
		echo "Setting permissions for file...\n";
		exec("chgrp \"domänen-benutzer\" {$to}");
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


$filename = basename($file, ".log");
if (! ends_with(".tar.gz", $file)) $filename .= ".tar.gz";

//check filename for forbidden symbols:
$matched = preg_match_all('/[^A-Za-z0-9\_\.]/', $filename);
if ($matched != 0)
{
	trigger_error("Filename \"{$filename}\" contains forbidden symbols only 'A-Z', 'a-z', '0-9', '_' and '.' are allowed.", E_USER_ERROR);
}


$file = "/mnt/SRV018/raw_data_archive/runs/{$filename}";
$to = "/mnt/storage1/raw_data/_restored/{$filename}";

if (file_exists($to)) trigger_error("The file $filename already exists in the outFolder!", E_USER_ERROR);
if (file_exists("/mnt/storage1/raw_data_archive/runs/{$filename}")) trigger_error("The file still exists in the archiving folder: /mnt/storage1/raw_data_archive/runs/{$filename}", E_USER_ERROR);

echo "Searching for $filename in archive.\n";

query($file);

$input = readline("Press ENTER, to start the restoration of the file. Press Strg+C to abort.");

restore($file, $to);


?>