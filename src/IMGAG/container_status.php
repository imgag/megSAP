<?php
/** 
	@page container_status
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("container_status", "Show status of containers.");
$parser->addOutfile("tsv", "Store output as TSV in addition to printing it.", true);
$parser->addFlag("check_md5", "Checks if MD5 checksums are correct (slow)");
extract($parser->parse($argv));

//init
$container_recipes_folder = repository_basedir()."/data/tools/container_recipes/";
$container_folder = get_path("container_folder");

//determine containers used in settings files
function add_containers($file, $display_name="")
{
	global $containers;
	
	$filename = repository_basedir()."/".$file;
	if (!file_exists($filename)) return;
	
	if ($display_name=="") $display_name = $file;
	
	foreach(file($filename) as $line)
	{
		$line = trim($line);
		if (!starts_with($line, "container_")) continue;
		if (starts_with($line, "container_folder")) continue;
		
		list($name, $version) = explode("=", $line, 2);
		$name = trim(strtr($name, ["container_"=>""]));
		$version = trim($version);
		$containers[$name][$version][] = $display_name;
	}
}
$containers = [];
add_containers("settings.ini");
add_containers("settings_nightly.ini", "nightly");
add_containers("settings.ini.default", "default");

//check containers
$used_files = [];
$filename = ($tsv=="") ? $parser->tempFile(".tsv") : $tsv;
$h = fopen2($filename, "w");
fputs($h, implode("\t", ["#tool", "tag", "used_in", "local_sif", "size_mb", "local_md5", "remote_sif", "remote_md5", "recipe", "log"])."\n");
foreach($containers as $tool => $tmp)
{
	foreach($tmp as $tag => $used_in)
	{
		$used_in = trim(implode(", ", $used_in));
		
		//check local SIF file
		$sif = $container_folder."/{$tool}_{$tag}.sif";
		$local_sif = file_exists($sif) ? "yes" : "no";
		$size_mb = $local_sif=="yes" ? number_format(filesize($sif)/1024/1024, 0, '.', '') : "";
		
		//check local MD5 file
		$md5 = $container_folder."/checksums/{$tool}_{$tag}.sif.md5";
		$local_md5 = file_exists($md5) ? "yes" : "no";
		$checksum = NULL;
		if ($check_md5 && $local_sif=="yes" && $local_md5=="yes")
		{
			list($stdout) = exec2("md5sum {$sif}");
			$checksum = explode(" ", $stdout[0])[0];
			$checksum2 = trim(file_get_contents($md5));
			if ($checksum!=$checksum2) $local_md5 .=" (incorrect)";
		}
		
		//check remote SIF file
		$remote_sif = url_exists("https://megsap.de/download/container/{$tool}_{$tag}.sif") ? "yes" : "no";
		
		//check remote MD5 file
		$checksum2 = "";
		$remote_md5 = url_exists("https://megsap.de/download/container/checksums/{$tool}_{$tag}.sif.md5", $checksum2) ? "yes" : "no";
		if ($check_md5 && !is_null($checksum) && $remote_md5=="yes")
		{
			if ($checksum!=$checksum2) $remote_md5 .=" (incorrect)";
		}
		
		//check recipe
		$tag_nodate = $tag;
		if (preg_match('/-20\d{6}$/', $tag)) $tag_nodate = substr($tag, 0, strlen($tag)-9); //remove date
		$recipe = "no";
		if (file_exists($container_recipes_folder."{$tool}_{$tag_nodate}.def"))
		{
			$recipe = "yes";
			$used_files["{$tool}_{$tag_nodate}.def"] = true;
		}
		if ($recipe=="no" && file_exists($container_recipes_folder."{$tool}_{$tag_nodate}_pull_command.txt"))
		{
			$recipe = "yes (docker)";
			$used_files["{$tool}_{$tag_nodate}_pull_command.txt"] = true;
		}
		
		//check log
		$log = "no";
		if (file_exists($container_recipes_folder."{$tool}_{$tag}.log.gz"))
		{
			$log = "yes";
			$used_files["{$tool}_{$tag}.log.gz"] = true;
		}
		
		fputs($h, implode("\t", [$tool, $tag, $used_in, $local_sif, $size_mb, $local_md5, $remote_sif, $remote_md5, $recipe, $log])."\n");
	}
}
fclose($h);

//display TSV in console
passthru("TsvTo -in {$filename} -format txt");
print "\n";

//check for extra recipes/logs that are not longer used
list($files) = exec2("ls {$container_recipes_folder}");
foreach($files as $file)
{
	if ($file=="backups") continue;
	if (starts_with($file, "ngs-bits_")) continue; //we never use release version during development...
	
	if (!isset($used_files[$file])) print "Note: Extra file: ".realpath($container_recipes_folder."/".$file)."\n";
}

//check if tools are used in pipeline
foreach($containers as $tool => $tmp)
{
	list($hits) = exec2("find src/ -name '*.php' | xargs egrep 'execApptainer.*\"{$tool}\"'", false);
	if (trim(implode("", $hits))=="")
	{
		print "Note: Unused tool: {$tool}\n";
	}
}

?>