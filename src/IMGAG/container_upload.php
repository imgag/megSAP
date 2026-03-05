<?php
/** 
	@page container_upload
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("container_upload", "Uploads an apptainer container to megsap.de and deploys it to the container repository.");
$parser->addString("tool", "Tool name", false);
$parser->addString("tag", "Tool tag or version.", true, "master");
$parser->addString("pw", "Password for megsap.de", true, "");
$parser->addFlag("no_upload", "Only deploys the container to the container repo without uploading it to megsap.de");
$parser->addFlag("copy", "Copy container instead of moving (for testing pipelines).");
extract($parser->parse($argv));

function deploy_container($log, $sif, $md5, $container_repo, $copy)
{
	if ($copy)
	{
		$desc = "Copying";
		$cmd = "cp";
	}
	else
	{
		$desc = "Moving";
		$cmd = "mv";
	}
	
	$recipe_folder = repository_basedir()."/data/tools/container_recipes/";
	print "{$desc} {$log} to {$recipe_folder}\n";
	exec2("{$cmd} {$log} {$recipe_folder}");
	if (!$copy) exec2("gzip -9 {$recipe_folder}/{$log}");
	print "{$desc} {$sif} to {$container_repo}/{$sif}\n";
	exec2("{$cmd} {$sif} {$container_repo}/{$sif}");
	print "{$desc} MD5 sum for $sif to {$container_repo}/checksums/{$sif}.md5\n";
	exec2("{$cmd} {$md5} {$container_repo}/checksums/{$md5}");
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
		trigger_error("Invalid input. Aborting!", E_USER_ERROR);
	}
}

//prevent undersore in tag
if (contains($tag, "_"))
{
	$tag = strtr($tag, "_", "-");
	print "Note: tool version/tag in container names must not contain underscore - using '{$tag}'\n";
}

//not ngs-bits master
if ($tool=="ngs-bits" && $tag=="master")
{
	print "FAILED: ngs-bits container built from current master is automatically deployed after building and should not be uploaded to megsap.de.\nAborting!\n";
	exit(1);
}

//get file names
exec("ls | grep '^" . preg_quote($tool."_".$tag, '/') . ".*\\.sif$'", $sif);
@$sif = trim($sif[0]);
if ($sif=="") trigger_error("SIF file not found!", E_USER_ERROR);
print "SIF: {$sif}\n";
exec("ls | grep '^" . preg_quote($tool."_".$tag, '/') . ".*\\.log$'", $log);
@$log = trim($log[0]);
if ($log=="") trigger_error("LOG file not found!", E_USER_ERROR);
print "LOG: {$log}\n";
exec("ls | grep '^" . preg_quote($tool."_".$tag, '/') . ".*\\.md5$'", $md5);
@$md5 = trim($md5[0]);
if ($md5=="") trigger_error("MD5 file not found!", E_USER_ERROR);
print "MD5: {$md5}\n";

//upload
if (!$no_upload)
{
	print "Did you test the container before uploading it to megsap.de?\n";
	print "(y/n)\n";

	$input = trim(fgets(STDIN));

	if (strtolower($input) === "y") 
	{
		if ($pw == "") trigger_error("Password for megsap.de must be set. Use parameter '-pw'", E_USER_ERROR);
		
		$user = trim(exec('whoami'));
		$web_container_dir = "{$user}@megsap.de:/mnt/nginx/public_html/download/container/";
		print "Uploading $sif and $md5 to megsap.de ...\n";

		//Check if container and checksum file are already present on megsap.de
		if (url_exists("https://megsap.de/download/container/{$sif}"))
		{
			$remote_md5 = "";
			if (!url_exists("https://megsap.de/download/container/checksums/{$md5}", $remote_md5))
			{
				//Ask if container should be uploaded if container is already on megsap.de but md5 sum is missing
				$message = "MD5 sum file missing for existing container {$sif} on megsap.de. Do you want to replace the existing container and upload the new MD5 sum file?\n";
				replace_remote_container($message, $sif, $pw, $md5, $web_container_dir);
			}
			else
			{
				//Check diff of local and remote md5
				if ($remote_md5==trim(file_get_contents($md5)))
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
			exec2("sshpass -p $pw scp -o PreferredAuthentications=password $sif {$web_container_dir}{$sif}");
			exec2("sshpass -p $pw scp -o PreferredAuthentications=password $md5 {$web_container_dir}checksums/{$md5}");
		}
	}
	elseif (strtolower($input) === "n") 
	{
		print "Aborting container upload...\n";
    	exit(1);
	} 
	else
	{
		trigger_error("Invalid input. Aborting!", E_USER_ERROR);
	}
}

//move/copy to container repository
$container_repo = "/mnt/storage2/megSAP/tools/apptainer_container/";
print "Deploying container {$sif} to {$container_repo}\n";
if (file_exists("{$container_repo}/{$sif}"))
{
	if (!file_exists("{$container_repo}/checksums/{$md5}"))
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
			deploy_container($log, $sif, $md5, $container_repo, $copy);
		}
		elseif (strtolower($input) === "n")
		{
			exit(1);
		} 
		else
		{
			trigger_error("Invalid input. Aborting!", E_USER_ERROR);
		}
	}
}
else 
{
	deploy_container($log, $sif, $md5, $container_repo, $copy);
}

?>