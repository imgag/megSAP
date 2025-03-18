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

//prevent undersore in tag
if (contains($tag, "_"))
{
	$tag = strtr($tag, "_", "-");
	print "Note: tool version/tag in container names must not contain underscore - using '{$tag}'\n";
}

$sif = "{$tool}_{$tag}.sif";
$log = "{$tool}_{$tag}.log";

if ($tool=="ngs-bits" && $tag=="master")
{
	print "FAILED: ngs-bits container built from current master is automatically deployed after building and should not be uploaded to megsap.de.\nAborting upload_apptainer_container.php.\n";
	exit(1);
}

if (!is_file($sif) || !is_file($log))
{
	trigger_error("Either $sif or $log not found. Check that both the container and the log file created by build_apptainer_container.php are present", E_USER_ERROR);
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
		print "Uploading $sif to megsap.de ...\n";
		exec2("sshpass -p $pw scp $sif megsap@megsap.de:/public_html/download/container/{$sif}");
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
	print "$sif already exists in $container_repo. Do you want to replace the existing container?\n";
	print "(y/n)\n";

	$input = trim(fgets(STDIN));

	if (strtolower($input) === "y") 
	{
		print "Moving $log to ".repository_basedir()."/data/tools/container_recipes/\n";
		exec2("mv {$log} ".repository_basedir()."/data/tools/container_recipes/");
		print "Moving $sif to {$container_repo}/{$sif}\n";
		exec2("mv {$sif} {$container_repo}/{$sif}");
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
else 
{
	print "Moving $log to ".repository_basedir()."/data/tools/container_recipes/\n";
	exec2("mv {$log} ".repository_basedir()."/data/tools/container_recipes/");
	print "Moving $sif to {$container_repo}/{$sif}\n";
	exec2("mv {$sif} {$container_repo}/{$sif}");
}
print "Deploying finished\n";

?>