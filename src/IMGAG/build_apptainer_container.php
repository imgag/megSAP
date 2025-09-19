<?php
/** 
	@page build_apptainer_container
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("build_apptainer_container", "Builds an apptainer container for a specified tool.");
$parser->addString("tool", "Tool to build a container for. Use 'which' to see available tools and their tags.", true, "which");
$parser->addString("tag", "Tag to build.", true, "master");
extract($parser->parse($argv));

//prevent undersore in tag
if (contains($tag, "_"))
{
	$tag = strtr($tag, "_", "-");
	print "Note: tool version/tag in container names must not contain underscore - using '{$tag}'\n";
}

//init
$date = date("Y-m-d H:i:s");
$datestr = date("Ymd");
$sif = "{$tool}_{$tag}-{$datestr}.sif";
$recipe_folder = "data/tools/container_recipes/";
$def = "{$recipe_folder}{$tool}_{$tag}.def";
$log = "{$tool}_{$tag}-{$datestr}.log";
$pull_file = "{$recipe_folder}{$tool}_{$tag}_pull_command.txt";
$container_recipes = repository_basedir()."/data/tools/container_recipes/";

//get available containers to build
if ($tool == "which") 
{
    $files = glob("$container_recipes/*.{def,txt}", GLOB_BRACE);
    $tools = [];

    foreach ($files as $file) 
	{
        $base = basename($file);

        if (preg_match('/^(.*?)_(.*?)\.def$/', $base, $m)) $tools[$m[1]][] = $m[2];
        elseif (preg_match('/^(.*?)_(.*?)_pull_command\.txt$/', $base, $m)) $tools[$m[1]][] = $m[2];
    }

    //format output
    printf("%-30s %s\n", "tools", "tags");
	printf("%'-30s %s\n", str_repeat("-", 5), str_repeat("-", 4));

    foreach ($tools as $t => $tags) 
	{
        foreach (array_unique($tags) as $tag) 
		{
            printf("%-30s %s\n", $t, $tag);
        }
    }

    exit();
}

//build
if (is_file($def))
{
	print "Building container {$sif} - in case of error see {$log}\n";

	// write server name, user and date to logfile as header
	$server = php_uname('n');
	$user = getenv("USER") ?: getenv("LOGNAME");
	

	$log_header = "Built on: {$server}\nBuilt by: {$user}\nBuild date: {$date}\n\n";
	file_put_contents($log, $log_header);
	exec2("apptainer build {$sif} {$recipe_folder}/{$tool}_{$tag}.def >> $log 2>&1");
}
elseif (is_file($pull_file))
{
	print "No definition file available for $tool, but found a pull command to download container. \n";
	exec("cat $pull_file", $pull_command);
	print "Pulling container {$sif} - in case of error see {$log}\n";

	// write server name, user and date to logfile as header
	$server = php_uname('n');
	$user = getenv("USER") ?: getenv("LOGNAME");

	$log_header = "Built on: {$server}\nBuilt by: {$user}\nBuild date: {$date}\nPull command: ".$pull_command[0]."\n\n";
	file_put_contents($log, $log_header);

	exec2($pull_command[0]." >> $log 2>&1");
	exec2("mv {$tool}_{$tag}.sif $sif");
}
else
{
	exec("ls -1 $recipe_folder | grep '^{$tool}.*\\.def$'", $def_available);
	print "Recipe '$def' not found! The following recipes are available for your tool: \n".implode("\n", $def_available)."\n";
	exit();
}

if ($tool != "ngs-bits" || $tag != "master")
{
	print "Calculating MD5 sum for $sif\n";
	exec2("md5sum -b {$sif} | cut -d ' ' -f1 > {$sif}.md5");
}

exec2("chmod 777 {$sif}");
print "Building container finished.\n";

if ($tool =="ngs-bits")
{
	//determine version
	list($stdout) = exec2("apptainer exec {$sif} MappingQC --version");
	$version = trim(strtr(implode("", $stdout), ["MappingQC"=>""]));
	print "ngs-bits version determined from container: {$version}\n";

	//determine final container name
	if($tag=="master")
	{
		$sif2 = "ngs-bits_master-".strtr($version, "_", "-").".sif";
		$log2 = "ngs-bits_master-".strtr($version, "_", "-").".log";

		//Calculate MD5 sum and move container and .md5 to container repo
		$container_repo = "/mnt/storage2/megSAP/tools/apptainer_container/";
		print "Calculating MD5 sum for $sif2\n";
		exec2("md5sum -b {$sif} | cut -d ' ' -f1 > {$container_repo}/checksums/{$sif2}.md5");

		print "Deploying container {$sif2} to {$container_repo}\n";
		exec2("mv {$sif} {$container_repo}/{$sif2}");
		print "Deploying finished\n";

		//Find old logs
		$old_logs = glob($container_recipes . "ngs-bits_master-*.log");

		//Remove any old logs (in case multiple exist)
		foreach ($old_logs as $old_log) 
		{
			print "Removing old logfile $old_log\n";
			unlink($old_log);
		}

		//Move new log
		print "Deploying new logfile to {$container_recipes}{$log2}\n";
		exec2("mv {$log} {$container_recipes}{$log2}");
	}
}

?>