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
       @brief Calculates a robust variant of the mean/stdev of a data array

       The calulation is performed by replacing @p num extreme values with the next smaller/larger values.

       @return an array with mean and standard deviation
       @ingroup statistics
*/
function winsored_statistics($data, $num)
{
	$n = count($data);

	sort($data);

	for ($i=0; $i<$num; ++$i)
	{
		$data[$i] = $data[$num];
		$data[$n-1-$i] = $data[$n-1-$num];
	}

	$mean = mean($data);
	$stdev = stdev($data, $mean);

	return array($mean, $stdev);
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

///Returns an array with all elements of @p heystack that contain @p needle.
function array_containing($haystack, $needle)
{
	return array_values(array_filter($haystack, function($var) use ($needle){ return strpos($var, $needle) !== false;}));
}

/**
 * @brief generates backup folder within data folder for a given pattern
 * @param type $data_folder folder to backup ('/' and emtpy are not allowed)
 * @param type $backup_pattern comma seprated lists of file extensions (e.g. *.tsv, *.log)
 * @param type $skip_pattern regular expression for patterns to exclude from backup (e.g. current log file)
 * @return type
 */
function backup($data_folder, $backup_pattern, $skip_pattern = NULL)
{
    //remove slashes
    $data_folder = trim($data_folder, "/");
    
    //skip empty folders or root folder
    if (empty($data_folder) || $data_folder == "/" || $data_folder == "/mnt")
    {
        trigger_error('Data folder not allowed for backup ("'.$data_folder.'").', E_USER_ERROR);
    }

    $old_files = glob($data_folder."/{".$backup_pattern."}", GLOB_BRACE);
    
    //filter duplicate files from backup_pattern
    $old_files = array_unique($old_files);
    
   //skip backup if no old files or only skipped files available.
	$counter = 0;
    foreach($old_files as $old_file)
    {
		//do not move files that match $skip_pattern
		if(!empty($skip_pattern) && preg_match($skip_pattern, $old_file))
		{
			$counter++;
		}
    }   
	if((count($old_files)-$counter) <= 0)
	{
		return;
	}

   //create backup_folder within $data_folder
	$old_folder = $data_folder."/old_".date('Y_m_d')."_".random_string(4)."/";
    if (is_dir($old_folder) || !mkdir($old_folder))
    {
        trigger_error('Could not create backup folder: '.$old_folder, E_USER_ERROR);
    }

    //find all files to backup in $data_folder
    foreach($old_files as $old_file)
    {
		//do not move files that match $skip_pattern
		if(empty($skip_pattern) || !preg_match($skip_pattern, $old_file))
		{
			rename($old_file, $old_folder.basename($old_file));
		}
    }
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
	
//adds a # in the first line of a file
function addCommentCharInHeader($filename)
{
	$file = fopen2($filename,"r+");
	$old_contents = file_get_contents($filename);
	if(starts_with($old_contents[0],"#")) return;
	fwrite($file,"#");
	fwrite($file,$old_contents);
	fclose($file);
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
?>