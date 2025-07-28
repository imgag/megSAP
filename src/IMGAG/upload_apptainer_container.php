<?php
/** 
	@page upload_apptainer_container
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("upload_apptainer_container", "Uploads an apptainer container to megsap.de and deploys it to the container repository.");
$parser->addString("tool", "Tool name", false);
$parser->addString("tag", "Tool tag (usually version).", true, "master");
$parser->addString("pw", "Password for megsap.de", true, "");
$parser->addFlag("no_upload", "Only deploys the container to the container repo without uploading it to megsap.de");
extract($parser->parse($argv));

//Functions
function remote_file_exists($remote_path, $pw) 
{
    $sftp_cmd = "echo 'ls " . escapeshellarg($remote_path) . "' | " .
                "sshpass -p " . escapeshellarg($pw) . " " .
                "sftp -q -oBatchMode=no -oStrictHostKeyChecking=no megsap@megsap.de 2>&1";

    exec($sftp_cmd, $output, $code);

    foreach ($output as $line) 
	{
        if (str_contains($line, 'No such file') || str_contains($line, 'not found')) 
		{
            return false;
        }
    }

    return $code === 0;
}

function deploy_container($log, $sif, $md5, $container_repo)
{

	print "Moving $log to ".repository_basedir()."/data/tools/container_recipes/\n";
	exec2("mv {$log} ".repository_basedir()."/data/tools/container_recipes/");
	print "Moving $sif to {$container_repo}/{$sif}\n";
	exec2("mv {$sif} {$container_repo}/{$sif}");
	print "Moving MD5 sum for $sif to {$container_repo}/checksums/{$sif}.md5\n";
	exec2("mv {$md5} {$container_repo}/checksums/$md5");
}

function replace_remote_container($message, $sif, $pw, $md5, $web_container_dir)
{
	$web_md5_file = "{$web_container_dir}checksums/{$md5}";
	$web_sif_file = "{$web_container_dir}{$sif}";

	print $message;
	print "(y/n)\n";
	$input = trim(fgets(STDIN));

	if (strtolower($input) === "y") 
	{
		print "Uploading container $sif...\n";
		exec2("sshpass -p $pw scp $sif {$web_sif_file}");
		exec2("sshpass -p $pw scp $md5 {$web_md5_file}");
	}
	elseif (strtolower($input) === "n")
	{
		exit(1);
	} 
	else
	{
		trigger_error("Invalid input. Aborting upload_apptainer_container.php.", E_USER_ERROR);
	}
}

//MAIN
//prevent undersore in tag
if (contains($tag, "_"))
{
	$tag = strtr($tag, "_", "-");
	print "Note: tool version/tag in container names must not contain underscore - using '{$tag}'\n";
}

//get file names
exec("ls | grep '^" . preg_quote($tool."_".$tag, '/') . ".*\\.sif$'", $sif);
exec("ls | grep '^" . preg_quote($tool."_".$tag, '/') . ".*\\.log$'", $log);
exec("ls | grep '^" . preg_quote($tool."_".$tag, '/') . ".*\\.md5$'", $md5);
$sif = reset($sif);
$log = reset($log);
$md5 = reset($md5);

if ($tool=="ngs-bits" && $tag=="master")
{
	print "FAILED: ngs-bits container built from current master is automatically deployed after building and should not be uploaded to megsap.de.\nAborting upload_apptainer_container.php.\n";
	exit(1);
}

if (!is_file($sif) || !is_file($log) || !is_file($md5))
{
	trigger_error("Either sif, log or md5 file not found. Check that the container, the log file and the MD5 sum file created by build_apptainer_container.php are present", E_USER_ERROR);
}

if (!$no_upload)
{
	print "Did you test the container before uploading it to megsap.de?\n";
	print "(y/n)\n";

	$input = trim(fgets(STDIN));

	if (strtolower($input) === "y") 
	{
		if ($pw === "")
		{
			trigger_error("Password for megsap.de must be set. Use parameter '-pw'", E_USER_ERROR);
		}
		$web_container_dir = "megsap@megsap.de:/public_html/download/container/";
		print "Uploading $sif and $md5 to megsap.de ...\n";

		//Check if container and checksum file are already present on megsap.de
		$remote_md5 = temp_file(".md5");
		exec("sshpass -p '$pw' scp -q {$web_container_dir}checksums/{$md5} $remote_md5 2>/dev/null", $out, $code_md5);

		if (remote_file_exists("/public_html/download/container/{$sif}", $pw))
		{
			if ($code_md5 != 0)
			{
				//Ask if container should be uploaded if container is already on megsap.de but md5 sum is missing
				$message = "MD5 sum file missing for existing container {$sif} on megsap.de. Do you want to replace the existing container and upload the new MD5 sum file?\n";
				replace_remote_container($message, $sif, $pw, $md5, $web_container_dir);
			}
			else
			{
				//Check diff of local and remote md5
				exec("diff $remote_md5 $md5", $output, $code);
				if ($code==0)
				{
					print "Container {$sif} with matching MD5 sum already uploaded to megsap.de - skipping upload.\n";
				} 
				else 
				{
					$message = "A Container {$sif} with differing checksum was found on megsap.de. Do you want to replace the existing container?\n";
					replace_remote_container($message, $sif, $pw, $md5, $web_container_dir);
				}
			}
		}
		else
		{
			print "Uploading container $sif...\n";
			exec2("sshpass -p $pw scp $sif {$web_container_dir}{$sif}");
			exec2("sshpass -p $pw scp $md5 {$web_container_dir}checksums/{$md5}");
		}
	}
	elseif (strtolower($input) === "n") 
	{
		print "Aborting container upload...\n";
    	exit(1);
	} 
	else
	{
		trigger_error("Invalid input. Aborting upload_apptainer_container.php.", E_USER_ERROR);
	}
}
else
{
	print "Note: uploading $sif to megsap.de will be skipped\n";
}

//move to container repository
$container_repo = "/mnt/storage2/megSAP/tools/apptainer_container/";
print "Deploying container {$sif} to {$container_repo}\n";
if (is_file("{$container_repo}/{$sif}"))
{
	if (!is_file("{$container_repo}/checksums/{$md5}"))
	{
		print "MD5 checksum missing for {$container_repo}/{$sif}, creating it now...\n";
		exec2("md5sum -b {$container_repo}/{$sif} | cut -d ' ' -f1 > {$container_repo}/checksums/{$md5}");
	}

	exec("diff $md5 {$container_repo}/checksums/{$md5}", $output, $code);
	if ($code==0)
	{
		print "Container $sif with matching MD5 sum already deployed - skipping deployment.\n";
	}
	else
	{
		print "A Container {$sif} with differing checksum was found in {$container_repo}. Do you want to replace the existing container?\n";
		print "(y/n)\n";
		$input = trim(fgets(STDIN));

		if (strtolower($input) === "y") 
		{
			deploy_container($log, $sif, $md5, $container_repo);
		}
		elseif (strtolower($input) === "n")
		{
			exit(1);
		} 
		else
		{
			trigger_error("Invalid input. Aborting upload_apptainer_container.php.", E_USER_ERROR);
		}
	}
}
else 
{
	deploy_container($log, $sif, $md5, $container_repo);
}
print "Deploying finished\n";

?>