<?php
/** 
	@page container_build
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("container_build", "Builds an apptainer container for a specified tool.");
$parser->addString("tool", "Tool to build a container for. Use 'which' to see available tools and their tags.", true, "which");
$parser->addString("tag", "Tool tag or version.", true, "master");
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

//show available containers to build
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
if (file_exists($def))
{
	print "Building container {$sif} - in case of error see {$log}\n";

	// write server name, user and date to logfile as header
	$server = php_uname('n');
	$user = getenv("USER") ?: getenv("LOGNAME");
	

	$log_header = "Built on: {$server}\nBuilt by: {$user}\nBuild date: {$date}\n\n";
	file_put_contents($log, $log_header);
	exec2("apptainer build --force {$sif} {$recipe_folder}/{$tool}_{$tag}.def >> $log 2>&1");
}
else if (file_exists($pull_file))
{
	print "No apptainer recipe exists, using Docker pull command for {$tool}_{$tag}.\n";
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
	exit(1);
}
exec2("chmod 777 {$sif}");
print "Building container finished.\n";

//calculate checksum
print "Calculating MD5 sum for $sif\n";
exec2("md5sum -b {$sif} | cut -d ' ' -f1 > {$sif}.md5");

//rename files for ngs-bits master (use git version in name instead of date)
if ($tool=="ngs-bits" && $tag=="master")
{
	//determine version
	list($stdout) = exec2("singularity exec {$sif} MappingQC --version");
	$version = trim(strtr(implode("", $stdout), ["MappingQC"=>""]));
	$version = strtr($version, "_", "-");
	print "ngs-bits version determined from container: {$version}\n";
	$tag = "{$tag}-{$version}";
	
	//rename files
	$sif2 = "ngs-bits_{$tag}.sif";
	exec2("mv {$sif} {$sif2}");
	exec2("mv {$sif}.md5 {$sif2}.md5");
	$log2 = "ngs-bits_{$tag}.log";
	exec2("mv {$log} {$log2}");
}

//deployment infos
print "\n";
print "Note: If you want to test the container in megSAP, you can deploy it locally using:\n";
print "> php src/IMGAG/container_upload.php -tool {$tool} -tag {$tag} -no_upload -copy\n";
print "For use of the container outside UKT, you have to upload the contain later!\n";

?>