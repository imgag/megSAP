<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

/**
	@brief Caclulates the factorial of a number.
	@ingroup helpers
*/
function factorial($number)
{ 
    if ($number < 2)
	{ 
        return 1; 
    }
	else
	{ 
        return $number * factorial($number-1); 
    }
}

/**
	@brief Swaps two variables.
	@ingroup helpers
*/
function swap(&$x, &$y)
{
    $tmp=$x;
    $x=$y;
    $y=$tmp;
}

/**
	@brief Bounds a value by the interval [min, max]
	@ingroup helpers
*/
function bound($v, $min, $max)
{
	if ($v<$min) $v=$min;
	if ($v>$max) $v=$max;
	return $v;
}

/**
	@brief Removes the newline character(s) from the end of a string
	@ingroup helpers
*/
function nl_trim($string)
{
	return rtrim($string, "\n\r");
}

/**
	@brief Executes a command and returns stdout output, stderr output and exit code.
	If @p abort_on_error is 'true' and the exit code not '0', execution of the script is aborted.
	@ingroup helpers
*/
function exec2($command, $abort_on_error = true)
{
    //start processing
	$proc = proc_open("bash -c \"set -o pipefail && ".$command."\"", array(1 => array('pipe','w'), 2 => array('pipe','w')), $pipes);
	if ($proc===false) trigger_error("Could not start process with exec2 function!\nCommand: {$command}", E_USER_ERROR);
	
	//get stdout, stderr and exit code
    $stdout = stream_get_contents($pipes[1]);
    fclose($pipes[1]);
    $stderr = stream_get_contents($pipes[2]);
    fclose($pipes[2]);
    $exit = proc_close($proc);
	
	//abort if requested
	if ($abort_on_error && $exit!=0)
	{
		trigger_error("Error while executing command: '$command'\nCODE: $exit\nSTDOUT: ".$stdout."\nSTDERR: ".$stderr."\n", E_USER_ERROR);
	}
	
	//return output
	return array(explode("\n", nl_trim($stdout)), explode("\n", nl_trim($stderr)), $exit);
}

/**
	@brief Returns if two ranges overlap
	@ingroup helpers
*/
function range_overlap($from1, $to1, $from2, $to2)
{
	return ($from2<=$from1 && $to2>=$from1) || ($from2<=$to1 && $to2>=$to1) || ($from2>=$from1 && $to2<=$to1);
}

/**
	@brief Intersects two regions.
	
	@return Array with start and end of the new range. FALSE if the input ranges did not overlap.
	
	@ingroup helpers
*/
function range_intersect($from1, $to1, $from2, $to2)
{
	if (!range_overlap($from1, $to1, $from2, $to2))
	{
		return false;
	}
	
	return array(max($from1, $from2), min($to1, $to2));
}

/**
	@brief Calculates the mean of a data array
	@return The mean value
	@ingroup statistics
*/
function mean(&$data)
{
	return array_sum($data)/count($data);
}

/**
	@brief Calculates the standard deviation of a data array
	@return The standard deviation value
	@ingroup statistics
*/
function stdev(&$data, $mean=null)
{
	if($mean==null)
	{
		$mean = mean($data);
	}

	$stdev = 0.0;
	foreach ($data as $d)
	{
		$stdev += ($d-$mean)*($d-$mean);
	}
	
	return sqrt($stdev/count($data));
}

/**
	@brief Calculates the median of a data array
	@return The median value
	@ingroup statistics
*/
function median($data)
{
	sort($data);
	
	$n = count($data);
	if ($n%2==0)
	{
		return ($data[$n/2-1] + $data[$n/2])/2;	
	}
	else
	{
		return $data[($n-1)/2];
	}
}

/**
	@brief Calculates the Median absolute deviation (a robust estimator for the standard diviation).
	@return The MAD value (multiplied by 1.4826 because we assume the data is normal-distributed).
	@ingroup statistics
*/
function mad($data, $median = null)
{
	if (!isset($median))
	{
		$median = median($data);
	}
	
	$res_abs = array();
	
	foreach ($data as $value)
	{
		$res_abs[] = abs($value - $median);
	}
	
	return 1.4826 * median($res_abs);
}

/**
	@brief Calculates the Pearson correlation coefficient of two samples.
	@return The Pearson correlation coefficient.
	@ingroup statistics
*/
function correlation(&$x, &$y)
{
	$n = count($x);
	if ($n!=count($y))
	{
		trigger_error("Samples must have the same size!", E_USER_ERROR);
	}
	
	$x_mean = mean($x);
	$y_mean = mean($y);
	
	$sum = 0.0;
	for($i=0; $i<$n; ++$i)
	{
		$sum += ($x[$i]-$x_mean) * ($y[$i]-$y_mean);
	}
	
	return $sum / $n / stdev($x, $x_mean) / stdev($y, $y_mean);
}

/**
	@brief Returns if a string starts with a certain prefix.
	@ingroup helpers
*/
function starts_with($string, $prefix)
{
	return substr($string, 0, strlen($prefix)) == $prefix;
}

/**
	@brief Returns if a string contains a certain substring.
	@ingroup helpers
*/
function contains($haystack, $needle)
{
	if (gettype($haystack)!="string") $haystack = strval($haystack);
	if (gettype($needle)!="string") $needle = strval($needle);
	
	return strpos($haystack,$needle)!== false;
}

/**
	@brief Returns if a string ends with a certain suffix.
	@ingroup helpers
*/
function ends_with($string, $suffix)
{
	return substr($string, -strlen($suffix)) == $suffix;
}

/**
	@brief Returns a writable temporary file name.
	@ingroup fileio
*/
function temp_file($suffix = "", $prefix = "tmp", $directory = null)
{
	if (is_null($directory))
	{
		$user = trim(exec('whoami'));
		$directory = sys_get_temp_dir()."/megSAP_user_{$user}/";
		create_directory($directory);
	}
	
	$filename = tempnam($directory, $prefix);
	if ($suffix!="")
	{
		unlink($filename);
	}
	return $filename.$suffix;
}

/*
	@brief Returns a writable temporary path.
	@ingroup fileio
 */
function temp_folder($prefix='megSAP_', $mode=0700)
{
	$user = trim(exec('whoami'));
	$user_tmp_dir = sys_get_temp_dir()."/megSAP_user_{$user}/";
	create_directory($user_tmp_dir);
    
    do
    {
		$path = $user_tmp_dir.$prefix.random_string(6);
    }
	while(file_exists($path) || !mkdir($path, $mode));
    
    return $path;
}

//Returns if the given file/directory is in the tmp folder
function is_in_temp_folder($filename)
{
	if (file_exists($filename)) $filename = realpath($filename);
	$temp = realpath(sys_get_temp_dir())."/";
	return starts_with($filename, $temp);
}


/**
	@brief Returns if an XML string is wellformed.
	@ingroup xml
*/
function xml_is_wellformed($string, &$messages)
{
	libxml_use_internal_errors(true);
	
	$doc = simplexml_load_string($string);
	
	$messages = libxml_get_errors();
	libxml_clear_errors();
	foreach ($messages as $message)
	{
		if ($message->level==LIBXML_ERR_ERROR || $message->level==LIBXML_ERR_FATAL)
		{
				return false;
		}
	}
	
	return true;
}


/**
	@brief Returns if an XML string matches a given schema.
	@ingroup xml
*/
function xml_matches_schema($string, $schema, &$messages)
{
	libxml_use_internal_errors(true);
	
	$xdoc = new DomDocument;
	$xdoc->LoadXML($string);
	
	$result = $xdoc->schemaValidateSource($schema);
			
	$messages = libxml_get_errors();
	libxml_clear_errors();
	
	return $result;
}


/**
	@brief Prints XML error and warning messages array (LibXMLError objects) to stdout or to a given array.
	
	If the input lines are given, the error position is printed as well.
	
	@ingroup xml
*/
function xml_print_messages($messages, $lines = null, $to_stdout = true)
{
	$output = array();
	foreach ($messages as $message)
	{
		$line = "";
		switch ($message->level)
		{
			case LIBXML_ERR_WARNING:
				$line .= "Warning";
				break;
			 case LIBXML_ERR_ERROR:
				$line .= "Error";
				break;
			case LIBXML_ERR_FATAL:
				$line .= "Fatal error";
				break;
		}
		$line .= " (".$message->code."): ".nl_trim($message->message)." in line $message->line, column $message->column.";
		$output[] = $line;
		
		if (isset($lines))
		{
			$output[] = "> ".nl_trim($lines[$message->line-1]);
			$output[] = "> ".str_pad ("|", $message->column-1, " ", STR_PAD_LEFT);
		}
	}
	
	if ($to_stdout)
	{
		print(implode("\n", $output));
	}
	
	return $output;
}

/**
	@brief Creates a random string of length @p length from the @p characters
	@ingroup helpers
*/
function random_string($length, $characters = "ABCDEFGHIJKLMNOPQRSTUVWQYZabcdefghijklmnopqrstuvwqyz0123456789")
{
	//init
	$out = "";
	mt_srand((double)microtime()*1000000);
	
	//create string
	for ($i=0;$i<$length;$i++)
	{
		$out .= $characters[mt_rand(0, strlen($characters)-1)];
	}
	return $out;
}

/**
	@brief Returns an ISO timestamp of the current time with or without microseconds
	@ingroup helpers
*/
function get_timestamp($add_ms=true)
{
	date_default_timezone_set('Europe/Berlin');
	$time = microtime(true);
	$date = date("Y-m-d\TH:i:s", $time);
	if ($add_ms) $date .= ".".substr($time - floor($time), 2, 4);
	return $date;
}

/**
	@brief Returns the human-readable time of a duration in seconds.
	
	@ingroup helpers
*/
function time_readable($duration)
{
	$output = array();
		
	//hours
	if ($duration>=3600)
	{
		$output[] = intdiv($duration, 3600)."h";
		$duration %= 3600;
	}
	
	//minutes
	if (count($output)!=0 || $duration>=60)
	{
		$output[] = intdiv($duration, 60)."m";
		$duration %= 60;
	}
	
	//seconds
	if (count($output)!=0)
	{
		$output[] = number_format($duration, 0)."s";
	}
	else
	{
		$output[] = number_format($duration, 4)."s";		
	}
	
    return implode(" ", $output);
}

/*
	@brief Returns a human-readable traceback string.
*/
function traceback()
{
	$output = array();
	$output[] = "Traceback:";
	
	$tb = debug_backtrace();
	foreach($tb as $index => $data)
	{
		//skip this function and the error handler
		if ($index<=1) continue;
		
		if (isset($data["file"]) && isset($data["line"]))
		{
			$output[] = "  ".$data["function"]." called at ".$data["file"].":".$data["line"];
		}
		else
		{
			$output[] = "  ".$data["function"];
		}
	}
	
	return implode("\n", $output)."\n";
}

///Returns the base path of the repository
function repository_basedir()
{
	return realpath(dirname(__FILE__)."/../../")."/";
}

///Returns the revision of the repository.
function repository_revision($prefix_repos_name=false)
{
	$tag_file = repository_basedir()."/megSAP_tag.txt";
	if (file_exists($tag_file))
	{
		$tag = trim(file_get_contents($tag_file));
	}
	else
	{
		$output = array();
		exec("cd ".repository_basedir()." && git describe --tags", $output);
		$tag = trim($output[0]);
	}
	
	return ($prefix_repos_name ? "megSAP " : "").$tag;
}

/*
	@brief Prints a data structure in human-readable HTML format.
*/
function print_h($value)
{
	print strtr(print_r($value,true), array("\n"=>"<br>"," "=>"&nbsp;"));
}

/*
	@brief Returns a subset array of $values with the elements defined by the indices given as parameters after $value.
	
	@note Instead of several indices, an array with indices can be provided as second element.
*/
function array_subset($values /*, $index1, index 2, ...*/)
{
	$args = func_get_args();
	if (is_array($args[1]))
	{
		$columns  = $args[1];
	}
	else
	{
		$columns = array_slice($args, 1);
	}
	
	$output = array();
	foreach($columns as $c)
	{
		$output[] = $values[$c];
	}
	return $output;
}

/*
	@brief Loads a tab-separated file without newline characters, empty lines and comment lines.
*/
function load_tsv($filename, $sep = "\t", $comment = "#")
{
	$output = array();
	
	$file = file($filename, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
	foreach($file as $line)
	{
		//skip comment lines
		if ($line[0]==$comment) continue;
		
		$output[] = explode($sep, $line);
	}
	
	return $output;
}

/**
	@brief Returns the common prefix of two strings.
*/
function common_prefix($s1, $s2)
{
    $i = 0;
    $max_len = min(strlen($s1), strlen($s2));
    while($i<$max_len && $s1[$i]==$s2[$i])
    {
		++$i;
    }
	return substr($s1, 0, $i);
}

/**
	@brief Returns the common suffix of two strings.
*/
function common_suffix($s1, $s2)
{
    $i = 0;
	$len_s1 = strlen($s1);
	$len_s2 = strlen($s2);
    $max_len = min($len_s1, $len_s2);
    while($i<$max_len && $s1[$len_s1-$i-1]==$s2[$len_s2-$i-1])
    {
		++$i;
    }
	
	if ($i==0) return "";
	
	return substr($s1, -$i);
}

///checks if directory exists and creates it if not
function create_directory($path)
{
	if (!is_dir($path))
	{
		mkdir($path);
	}
}

///Own implementation of is_writable, checks if existing file is writable or if writable file can be created otherwise
function is_writable2($file)
{
	$existed = file_exists($file);

	$handle = @fopen($file, 'a');
	if ($handle===false)
	{
		$writable = false;
	}
	else
	{
		$writable = true;
		fclose($handle);
	}
	
	if (!$existed && file_exists($file))
	{
		unlink($file);
	}
	
	return $writable;
}

function sort_vcf_comments($comments_to_sort)
{
	//group comments
	$groups = array(
		"other"=>array(),
		"#contig"=>array(),
		"#INFO"=>array(),
		"#FILTER"=>array(),
		"#FORMAT"=>array(),
		"#ALT"=>array(),
		"#assembly"=>array(),
		"#SAMPLE"=>array()
		);
	foreach($comments_to_sort as $comment)
	{
		list($i) = explode("=", $comment, 2);
		if (starts_with($i, "##")) $i = substr($i, 1); //leading "##" => "#"
		if(isset($groups[$i]))
		{
			$groups[$i][] = $comment;
		}
		else
		{
			$groups["other"][] = $comment;
		}
	}
	
	//sort by group / text
	$sorted = array();
	foreach($groups as $group => $comments)
	{
		//"other" contains fileformat and similar headers, which must not change order!
		//"#SAMPLE" order must not change, otherwise the order of files shown in IGV is random
		if ($group!="other" && $group!="#SAMPLE") 
		{
			sort($comments, SORT_FLAG_CASE|SORT_STRING);
		}
		$sorted = array_merge($sorted, $comments);
	}

	return $sorted;
}

//open an file and returns the file handle. Throws an error if it fails!
function fopen2($filename, $mode)
{
	$handle = fopen($filename, $mode);
	
	if ($handle===false)
	{
		trigger_error("Could not open file '$filename' in mode '$mode'!", E_USER_ERROR);
	}
	
	return $handle;
}
	
//open a GZ file and returns the file handle. Throws an error if it fails!
function gzopen2($filename, $mode)
{
	$handle = gzopen($filename, $mode);
	
	if ($handle===false)
	{
		trigger_error("Could not open GZ file '$filename' in mode '$mode'!", E_USER_ERROR);
	}
	
	return $handle;
}

//functions to encode and decode VCF INFO values using URL encoding
$vcf_encode_mapping = array("%" => "%25", "\t" => "%09", "\n" => "%0A", "\r" => "%0D", " " => "%20", "," => "%2C", ";" => "%3B", "=" => "%3D", "&" => "%26");
$vcf_decode_mapping = array_reverse(array_flip($vcf_encode_mapping));

// encode string using URL encoding
function vcf_encode_url_string($input_string)
{
	//encode string
	return strtr($input_string, $GLOBALS['vcf_encode_mapping']);
}
// decode URL encoded string
function vcf_decode_url_string($encoded_string)
{
	// define decode mapping:
	// (also replaces newlines and tabs with spaces)
	$additional_mapping = array("\t" => " ", "\n" => " ", "\r" => " ");
	return strtr(strtr($encoded_string, $GLOBALS['vcf_decode_mapping']), $additional_mapping);
}

//returns the real file name of a symlink to a file, or the given file name if it is not a symlink.
function resolve_symlink($filename)
{
	//no file > return input
	if (!is_file($filename)) return $filename;
	
	//no link > return input
	if (!is_link($filename)) return $filename;
	
	return realpath($filename);
}
/**
	@brief Executes a command inside a given Apptainer container and returns an array with STDOUT, STDERR and exit code.
*/
function execApptainer($container, $command, $parameters, $in_files=[], $out_folders=[], $command_only=false, $return_for_toolbase=false, $abort_on_error=true, $gpu_container=false)
{
	//check input
	if (is_array($command))
	{
		trigger_error("Error in 'execApptainer': command cannot be array: ".implode("\n", $command), E_USER_ERROR);
	}
	if (is_array($parameters))
	{
		trigger_error("Error in 'execApptainer' of '{$command}': parameters cannot be array: ".implode("\n", $parameters), E_USER_ERROR);
	}
	if(contains($command." ".$parameters, "|"))
	{
		trigger_error("Error in 'execApptainer' of '{$command}': command must not contain pipe symbol '|': $command $parameters", E_USER_ERROR);
	}
	if (!is_array($in_files))
	{
		trigger_error("Error in 'execApptainer' of '{$command}': in_files must be array!", E_USER_ERROR);
	}
	if (!is_array($out_folders))
	{
		trigger_error("Error in 'execApptainer' of '{$command}': out_folders must be array!", E_USER_ERROR);
	}
	
	//apptainer arguments
	$apptainer_args = [];
	$apptainer_args[] = "--no-mount home,cwd";
	$apptainer_args[] = "--cleanenv";

	if ($container=="samtools" || $container=="methylartist" || $container=="straglrOn") //samtools (and methylartist, straglrOn) need REF_CACHE/REF_PATH to avoid re-initializing the reference cache inside the container every call: http://www.htslib.org/doc/samtools.html#ENVIRONMENT_VARIABLES, https://www.htslib.org/workflow/cram.html#The%20REF_PATH%20and%20REF_CACHE 
	{
		$ref_cache_folder = get_path("local_data")."/samtools_ref_cache/";
		if (!file_exists($ref_cache_folder)) $ref_cache_folder = get_path("data_folder")."/genomes/samtools_ref_cache/";
		$apptainer_args[] = "--env REF_CACHE={$ref_cache_folder}/%2s/%2s/%s";
		$apptainer_args[] = "--env REF_PATH={$ref_cache_folder}/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s";
	}

	//set working dir to /tmp (otherwise the user home is used as default WD. But it is not mounted and this generates a warning)
	$apptainer_args[] = "--pwd=/tmp/";
		
	//to run a gpu supported apptainer container you need the --nv flag
	if ($gpu_container)
	{
		$apptainer_args[] = "--nv";
		$apptainer_args[] = "--env TF_FORCE_GPU_ALLOW_GROWTH=true"; //to avoid OOM error when using container with GPU support
	}

	//if ngs-bits container is executed the settings.ini is mounted into the container during execution 
	$bind_paths = array();
	if($container=="ngs-bits")
	{
		$ngsbits_local = get_path("ngs-bits_local", false);
		if ($ngsbits_local != "")
		{
			trigger_error("Using local ngs-bits installation specified in settings.ini.", E_USER_NOTICE);
			$ngsbits_command = $ngsbits_local.$command." ".$parameters;
			if ($command_only)
			{
				return $ngsbits_command;
			}
			else
			{
				list($stdout, $stderr, $return) = exec2($ngsbits_command, $abort_on_error);
				return array($stdout, $stderr, $return);
			}
		}
		
		//ngs-bits settings file cannot be in repo if running in container, so we put it into the data folder then
		$ngsbits_settings_loc = get_path("megSAP_container_used") ? get_path("data_folder")."/tools/ngsbits_settings.ini" : repository_basedir()."/data/tools/ngsbits_settings.ini";
		
		//ngs-bits settings file missing > create it
		if (!file_exists($ngsbits_settings_loc) || (file_exists(repository_basedir()."/settings.ini") && filemtime($ngsbits_settings_loc)<filemtime(repository_basedir()."/settings.ini")))
		{
			trigger_error("ngs-bits settings file at '{$ngsbits_settings_loc}' is missing or outdated. Updating it...", E_USER_NOTICE);
			$output = [];
			
			//reference genome
			$output[] = "reference_genome = ".genome_fasta("GRCh38");
			
			//NGSD credentials (production instance)
			$output[] = 'ngsd_host = "'.get_db('NGSD', 'db_host', ''). '"';
			$output[] = "ngsd_port = ".get_db('NGSD', 'db_port', '');
			$output[] = 'ngsd_name = "'.get_db('NGSD', 'db_name', ''). '"';
			$output[] = 'ngsd_user = "'.get_db('NGSD', 'db_user', ''). '"';
			$output[] = 'ngsd_pass = "' . get_db('NGSD', 'db_pass', '') . '"';
			
			//NGSD credentials (test instance)
			$output[] = 'ngsd_test_host = "'.get_db('NGSD_TEST', 'db_host', ''). '"';
			$output[] = "ngsd_test_port = ".get_db('NGSD_TEST', 'db_port', '');
			$output[] = 'ngsd_test_name = "'.get_db('NGSD_TEST', 'db_name', ''). '"';
			$output[] = 'ngsd_test_user = "'.get_db('NGSD_TEST', 'db_user', ''). '"';
			$output[] = 'ngsd_test_pass = "'.get_db('NGSD_TEST', 'db_pass', ''). '"';
			
			//GenLab credentials
			$output[] = "genlab_mssql = ".(get_path("genlab_mssql") ? "true" : "false");
			$output[] = 'genlab_host = "'.get_path("genlab_host"). '"';
			$output[] = "genlab_port = ".get_path("genlab_port");
			$output[] = 'genlab_name = "'.get_path("genlab_name"). '"';
			$output[] = 'genlab_user = "'.get_path("genlab_user"). '"';
			$output[] = 'genlab_pass = "'.get_path("genlab_pass"). '"';


			//project folders
			$project_folder = get_path("project_folder", false);
			$output[] = "projects_folder_diagnostic = ".(is_array($project_folder) ? $project_folder['diagnostic'] : $project_folder."/diagnostic/");
			$output[] = "projects_folder_research = ".(is_array($project_folder) ? $project_folder['research'] : $project_folder."/research/");
			$output[] = "projects_folder_test = ".(is_array($project_folder) ? $project_folder['test'] : $project_folder."/test/");
			$output[] = "projects_folder_external = ".(is_array($project_folder) ? $project_folder['external'] : $project_folder."/external/");
			
			//data folder
			$output[] = "data_folder = ".get_path("data_folder", false);
									
			$written = file_put_contents($ngsbits_settings_loc, implode("\n", $output));
			if($written===false) trigger_error("Could not write ngs-bits settings file: $ngsbits_settings_loc", E_USER_ERROR);

			//try to change permissions
			exec2("chmod 777 {$ngsbits_settings_loc}", false);

		}
		
		$parameters = "--settings {$ngsbits_settings_loc} ".$parameters;
		$bind_paths[] = $ngsbits_settings_loc.":".$ngsbits_settings_loc;
	}
	
	//get container (preferably from local folder)
	$container_version = get_path("container_{$container}");
	$container_file = "{$container}_{$container_version}.sif";
	$local_container = get_path("local_data")."/container/{$container_file}";
	$network_container = get_path("container_folder")."/{$container_file}";

	if (file_exists($local_container))
	{
		$container_path = $local_container;
	}
	else if (file_exists($network_container))
	{
		$container_path = $network_container;
	}
	else
	{
		trigger_error("Apptainer container '{$container_file}' neither found in '{$local_container}' nor in '{$network_container}'", E_USER_ERROR);
	}
	
	//determine bind paths from input and output files
	if (!get_path("megSAP_container_used"))
	{
		foreach($in_files as $file)
		{
			if (is_dir($file)) 
			{
				$filepath = realpath($file);
			} 
			else 
			{
				$filepath = dirname(realpath($file));
			}
			if(!in_array($filepath.":".$filepath, $bind_paths)) $bind_paths[] = $filepath.":".$filepath; 
		}

		foreach($out_folders as $folder)
		{
			//check it is a folder
			if (is_file($folder)) trigger_error("{$container}: Only folders can be bound as output parameters. '{$folder}' is a file!", E_USER_ERROR);
			
			// ! Bind parent folder of given output in case the out-folder doesn't exist yet. (Folder can't be bound if it doesn't exist)
			$filepath = realpath(dirname($folder));

			if(!in_array($filepath.":".$filepath, $bind_paths)) $bind_paths[] = $filepath.":".$filepath; 
		}
	}

	//handle binding of temp folder for SigProfilerExtractor
	if($container === "SigProfilerExtractor")
	{
		$templates_dir = temp_folder();
		$cosmic_templates_dir = temp_folder();
		$templates_bind_path = "{$templates_dir}:/opt/Python3.8.10/lib/python3.8/site-packages/sigProfilerPlotting/templates/";
		$cosmic_templates_bind_path = "{$cosmic_templates_dir}:/opt/Python3.8.10/lib/python3.8/site-packages/SigProfilerAssignment/DecompositionPlots/CosmicTemplates";

		if (get_path("megSAP_container_used"))
		{
			$apptainer_args[] = "-B {$templates_bind_path},{$cosmic_templates_bind_path}";
		}
		else
		{
			$bind_paths[] = $templates_bind_path;
			$bind_paths[] = $cosmic_templates_bind_path;
		}
	}
	
	//clean up empty paths
	for ($i=0; $i<count($bind_paths); ++$i)
	{
		if ($bind_paths[$i]==":")
		{
			unset($bind_paths[$i]);
		}
	}
	
	//check bind paths
	if(!empty($bind_paths) && !get_path("megSAP_container_used"))
	{
		foreach($bind_paths as $path)
		{
			//remove optional path option
			$path = explode(":", $path)[0];
			if(!file_exists($path))
			{
				trigger_error("Bind path '{$path}' does not exist!", E_USER_ERROR);
			}
		}
		$apptainer_args[] = "-B ".implode(",", $bind_paths);
	}
	else if (get_path("megSAP_container_used"))
	{
		$apptainer_args[] = "-B /megSAP/";
	}

	//compose Apptainer command
	$apptainer_command = "singularity exec ".implode(" ", $apptainer_args)." {$container_path} {$command} {$parameters}";

	if ($container=="deepvariant")
	{
		$apptainer_command = "TF_CPP_MIN_LOG_LEVEL=2 $apptainer_command";
	}
	
	//if command only option is true, only the apptainer command is being returned, without execution
	if($command_only) 
	{
		$apptainer_command = "singularity -q exec ".implode(" ", $apptainer_args)." {$container_path} {$command} {$parameters}";
		return $apptainer_command;
	}

	//return apptainer command, bind paths and the container path for the execApptainer function in ToolBase.php
	if($return_for_toolbase)
	{
		return array($apptainer_command, $bind_paths, $container_path, $container_version);
	}

	//execute call - pipe stdout/stderr to file
	$proc = proc_open($apptainer_command, array(1 => array('pipe','w'), 2 => array('pipe','w')), $pipes);
	if ($proc===false) trigger_error("Could not start process with ToolBase::execApptainer function!\nContainer: {$container_path}\nCommand: {$command}\nParameters: {$parameters}", E_USER_ERROR);
	
	//get stdout, stderr and exit code
	$stdout = stream_get_contents($pipes[1]);
	fclose($pipes[1]);
	$stderr = stream_get_contents($pipes[2]);
	fclose($pipes[2]);
	$exit = proc_close($proc);
	
	//abort if requested
	if ($abort_on_error && $exit!=0)
	{
		trigger_error("Error while executing command: '$apptainer_command'\nCODE: $exit\nSTDOUT: ".$stdout."\nSTDERR: ".$stderr."\n", E_USER_ERROR);
	}
	
	// Clean stderr by removing specific warning lines
	$stderr_lines = explode("\n", nl_trim($stderr));
	$filtered_stderr_lines = array_filter($stderr_lines, function($line) 
	{
		return stripos($line, "WARNING: Error changing the container working directory") === false;
	});

	// Return cleaned output
	return array(explode("\n", nl_trim($stdout)), array_values($filtered_stderr_lines), $exit);
}

// convert chromosome to NC number or the other way round.
function chr2NC($key, $revert=false, $throw_if_not_found=false)
{
	// init mapping tables (once)
	static $chr2nc = [
		"chr1" => "NC_000001.11",
		"chr2" => "NC_000002.12",
		"chr3" => "NC_000003.12",
		"chr4" => "NC_000004.12",
		"chr5" => "NC_000005.10",
		"chr6" => "NC_000006.12",
		"chr7" => "NC_000007.14",
		"chr8" => "NC_000008.11",
		"chr9" => "NC_000009.12",
		"chr10" => "NC_000010.11",
		"chr11" => "NC_000011.10",
		"chr12" => "NC_000012.12",
		"chr13" => "NC_000013.11",
		"chr14" => "NC_000014.9",
		"chr15" => "NC_000015.10",
		"chr16" => "NC_000016.10",
		"chr17" => "NC_000017.11",
		"chr18" => "NC_000018.10",
		"chr19" => "NC_000019.10",
		"chr20" => "NC_000020.11",
		"chr21" => "NC_000021.9",
		"chr22" => "NC_000022.11",
		"chrX" => "NC_000023.11",
		"chrY" => "NC_000024.10",
		"chrMT" => "NC_012920.1"
	];
	$nc2chr = array_flip($chr2nc);
	
	if (!$revert && isset($chr2nc[$key])) return $chr2nc[$key];
	if ($revert && isset($nc2chr[$key])) return $nc2chr[$key];
	
	if ($throw_if_not_found) trigger_error(__FUNCTION__.": cannot convert '{$key}'!", E_USER_ERROR);
	return "";
}

//returns which container platform is used: 'apptainer' or 'singularity'
function container_platform($add_version=false)
{
	list($stdout) = exec2("singularity --version");
	$stdout = strtolower(trim(implode(' ', $stdout)));
	$stdout = strtr($stdout, ["version "=>""]);
	
	$parts = explode(' ', $stdout);
	$platform = $parts[0];
	if ($platform!="apptainer" && $platform!="singularity") trigger_error("Invalid container platform: {$platform}", E_USER_ERROR);
	
	if ($add_version)
	{
		$platform .= ' '.$parts[1];
	}
	
	return $platform;
}

//check if a URL exists. If $content is set, it returns the content of the file.
function url_exists($url, &$content = null)
{
	//options
    $options = [
		CURLOPT_RETURNTRANSFER  => true,
        CURLOPT_FOLLOWLOCATION  => true,
        CURLOPT_TIMEOUT         => 5, //5s
        CURLOPT_CONNECTTIMEOUT  => 5, //5s
        CURLOPT_SSL_VERIFYPEER  => true,
        CURLOPT_SSL_VERIFYHOST  => 2,
    ];
	if (is_null($content)) $options[CURLOPT_NOBODY] = true;
	
	//exec
    $ch = curl_init($url);
	curl_setopt_array($ch, $options);
    $result = curl_exec($ch);
	
	//handle error
    if (curl_errno($ch))
	{
        curl_close($ch);
		if (!is_null($content)) $content = "";
        return false;
    }
	
	//get http code
    $http_code = curl_getinfo($ch, CURLINFO_HTTP_CODE);
    curl_close($ch);

    if ($http_code < 200 || $http_code >= 400)
	{
		if (!is_null($content)) $content = "";
        return false;
    }
	
	if (!is_null($content)) $content = trim($result);
	
	return true;
}
?>
