<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// Parse command line arguments
$parser = new ToolBase("container_setup", "Create local copies of frequently used Apptainer containers to prevent delays caused by accessing them from a network directory.");

extract($parser->parse($argv));

$network_folder = get_path("container_folder");
$local_folder = get_path("local_container_folder");
$rsync = "rsync --size-only --recursive --no-perms --no-acls --omit-dir-times --no-group --no-owner --chmod=ugo=rwX --copy-links";

print "### Copy apptainer containers ###\n";
print "from: {$network_folder}\n";
print "to  : {$local_folder}\n";
print "\n";

// Check if the network container folder exists
if (!file_exists($network_folder))
{
    trigger_error("Container folder not found: {$network_folder}. The Apptainer containers may not have been downloaded yet.", E_USER_ERROR);
}

// Create local container folder
if (!file_exists($local_folder))
{
	if (!mkdir($local_folder))
	{
		trigger_error("Could not create local data folder '{$local_folder}'!", E_USER_ERROR);
	}
	if (!chmod($local_folder, 0777))
	{
		trigger_error("Could not change privileges of local data folder '{$local_folder}'!", E_USER_ERROR);
	}
}

print "Copying apptainer containers...\n";
// Get list of apptainer containers from network directory
list($container_files) = exec2("ls {$network_folder}/*");

foreach ($container_files as $container_file) 
{
    $base = basename($container_file);
    $local_container_path = $local_folder . $base;

    // Split filename into toolname and version using underscore as the delimiter
    $parts = explode('_', $base);
    
    if (count($parts) >= 2) 
	{
        $toolname = $parts[0];
        
        $toolversion_with_extension = $parts[1];
        $toolversion = str_replace('.sif', '', $toolversion_with_extension);

        print "  checking container '$toolname' version '$toolversion'.\n";

        // Check if a local container with the same toolname exists
        $local_container_pattern = $local_folder . "{$toolname}_*.sif";
        $local_container_files = glob($local_container_pattern);

        if (!empty($local_container_files)) 
		{
            $local_base = basename($local_container_files[0]);
            $local_parts = explode('_', $local_base);
            
            if (count($local_parts) >= 2) 
			{
                $local_toolname = $local_parts[0];
                $local_toolversion = str_replace('.sif', '', $local_parts[1]);

                // Compare versions
                if ($toolversion !== $local_toolversion) 
				{
                    print "    newer version found: '$base'. Replacing local version '$local_base'.\n";
                    // Remove the old version
                    unlink($local_container_files[0]);
                    // Copy the new version
                    list($stdout, $stderr) = exec2("{$rsync} {$network_folder}{$base} {$local_folder}{$base}");
                } 
				else 
				{
                    print "    local version '$local_toolversion' is up to date for '$toolname'.\n";
                    continue;
                }
            }
        }
		else 
		{
            print "    copying new container '$base'.\n";
            list($stdout, $stderr) = exec2("{$rsync} {$network_folder}{$base} {$local_folder}{$base}");
        }

        foreach (array_merge($stdout, $stderr) as $line) 
		{
            $line = trim($line);
            if ($line == "") continue;
            print "    $line\n";
        }

        // Set permissions on the new local copy
        if (!chmod($local_container_path, 0777)) 
		{
            trigger_error("Could not change privileges of local folder '{$local_container_path}'!", E_USER_ERROR);
        }
    }
	else
	{
        print "    filename '$base' does not match expected pattern. Skipping...\n";
    }
}

?>