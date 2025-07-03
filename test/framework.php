<?php
require_once(dirname(__FILE__)."/../src/Common/all.php");

/// Returns the tool test data folder
function data_folder()
{
	return repository_basedir()."/test/data/";
}

/// Returns the tool test data folder for DB tests
function data_db_folder()
{
	return repository_basedir()."/test/data_db/";
}

/// Returns the tool test output folder
function output_folder()
{
	$folder = repository_basedir()."/test/data_out/".basename($_SERVER['SCRIPT_FILENAME'], ".php")."/";
	if (!file_exists($folder)) mkdir($folder);
	return $folder;
}

/// Returns the source folder
function src_folder()
{
	return repository_basedir()."/src/";
}

/// Converts a variable to a string (special handing of arrays)
function readable($data)
{
	if (!is_array($data))
	{
		if (is_bool($data))
		{
			return $data? "true": "false";
		}
		
		return $data;
	}
	
	return "array(".implode(",", array_map("readable", $data)).")";
}

/// Deletes all files from the output folder
function clear_output_folder()
{
	$dir = output_folder();
	if (file_exists($dir))
	{
		exec2("rm -rf $dir/*");
	}
}

///  Starts a test
function start_test($name)
{
	clear_output_folder();
	
	//checking for debug mode (argument '-d' on command line)
	$GLOBALS["debug"] = false;
	$args = getopt("d");
	if (isset($args["d"])) $GLOBALS["debug"] = true;
	$args = getopt("v");
	if (isset($args["v"])) $GLOBALS["debug"] = true;
	
	print "START $name\n";
	$GLOBALS["failed"] = 0;
	$GLOBALS["passed"] = 0;
	$GLOBALS["test_name"] = $name;
	$GLOBALS["start_time"] = microtime(true);
}

function check_test_started()
{
	if (!isset($GLOBALS["test_name"]))
	{
		trigger_error("Invalid test ".$_SERVER['SCRIPT_FILENAME'].": 'start_test' function not called!", E_USER_ERROR);
	}
}

/// Performs a check within a test
function check($observed, $expected, $delta = null)
{
	check_test_started();
	
	if (isset($delta))
	{
		$delta_string = ", delta='$delta'";
		if (is_array($observed))
		{
			$passed = true;
			for ($i = 0; $i < count($observed); ++$i)
			{
				$passed = $passed && (abs($observed[$i] - $expected[$i]) <= $delta);
			}
		}
		else
		{
			$passed = abs($observed - $expected) <= $delta;
		}
	}
	else
	{
		$delta_string = "";
		$passed = $observed==$expected;
	}
	
	if (!$passed)
	{
		$result = "FAILED (expected '".readable($expected)."', got '".readable($observed)."'$delta_string)";
		++$GLOBALS["failed"];
	}
	else
	{
		$result = "PASSED";
		++$GLOBALS["passed"];
	}
	
	$bt = debug_backtrace();
	$caller = array_shift($bt);
	$file = basename($caller["file"]);
	$line = $caller["line"];
	print "  - $file:$line $result\n";
}

///Removes lines that contain any string from the @p ignore_strings array 
function remove_lines_containing($filename, $ignore_strings)
{
	//file missing > nothing to do...
	if (!file_exists($filename)) return;

	//make sure ignore strings is an array
	if (!is_array($ignore_strings)) $ignore_strings = [$ignore_strings];
	
	//load input
	$file = file($filename);
	
	$h = fopen2($filename, "w");
	for($i=0; $i<count($file); ++$i)
	{
		$ignore = false;
		foreach($ignore_strings as $needle)
		{
			if (contains($file[$i], $needle))
			{
				$ignore = true;
				break;
			}
		}
		if (!$ignore)
		{
			fwrite($h, $file[$i]);
		}
	}
	fclose($h);
}

///Checks whether a file exists
function check_file_exists($in_file)
{
	check_test_started();

	if(file_exists($in_file))
	{
		$result = "PASSED";
		++$GLOBALS["passed"];
	}
	else
	{
		$result = "FAILED";
		++$GLOBALS["failed"];
	}
	
	$bt = debug_backtrace();
	$caller = array_shift($bt);
	$file = basename($caller["file"]);
	$line = $caller["line"];
	print "  - $file:$line $result\n";
}


///Checks whether a expression evaluates to true
function check_true($expression)
{
	check_test_started();

	if($expression)
	{
		$result = "PASSED";
		++$GLOBALS["passed"];
	}
	else
	{
		$result = "FAILED";
		++$GLOBALS["failed"];
	}
	
	$bt = debug_backtrace();
	$caller = array_shift($bt);
	$file = basename($caller["file"]);
	$line = $caller["line"];
	print "  - $file:$line $result\n";
}

///Checks whether a string contains a substring
function check_contains($haystack, $needle)
{
	check_test_started();

	if(contains($haystack, $needle))
	{
		$result = "PASSED";
		++$GLOBALS["passed"];
	}
	else
	{
		$result = "FAILED";
		++$GLOBALS["failed"];
		if ($GLOBALS["debug"]) print "    '{$needle}' not contained in '{$haystack}'\n";
	}
	
	$bt = debug_backtrace();
	$caller = array_shift($bt);
	$file = basename($caller["file"]);
	$line = $caller["line"];
	print "  - $file:$line $result\n";
}

//checks whether column with a certain name (that is: entry in 1st row) exists in a TSV file
function check_column_exists($filename,$col_names)
{
	check_test_started();

	$logfile = $filename ."_diff";
	$report_errors = "";
	
	$headers = file($filename)[0];
	$headers = str_replace("#","",$headers);
	$headers  = explode("\t",$headers);
	
	$passed = true;
	foreach($col_names as $col_name)
	{
		if(!in_array($col_name,$headers))
		{
			$report_errors .= "Column ".$col_name." does not exist in file ".$filename.".\n";
			$passed = false;
		}
	}

	if($passed)
	{
		$result = "PASSED";
		++$GLOBALS["passed"];
	}
	else
	{
		file_put_contents($logfile,$report_errors);
		$result = "FAILED";
		++$GLOBALS["failed"];
	}	
	
	$bt = debug_backtrace();
	$caller = array_shift($bt);
	$file = basename($caller["file"]);
	$line = $caller["line"];
	print "  - $file:$line $result\n";
	
}

/// Performs an equality check on files. Optionally, header lines starting with '#' can be compared as well.
function check_file($out_file, $reference_file, $compare_header_lines = false)
{
	check_test_started();
	
	$logfile = $out_file."_diff";
	
	//zdiff
	if (ends_with($out_file, ".gz") && ends_with($reference_file, ".gz"))
	{
		$extras = "";
		if (!$compare_header_lines) $extras .= " -I^[#@]";
		exec("zdiff $extras -b $reference_file $out_file > $logfile 2>&1", $output, $return);
		$passed = ($return==0 && count(file($logfile))==0);
	}
	//zipcmp
	elseif (ends_with($out_file, ".zip") && ends_with($reference_file, ".zip"))
	{
		exec("zipcmp $reference_file $out_file > $logfile 2>&1", $output, $return);
		$passed = ($return==0 && count(file($logfile))==0);
	}
	//bam (diff 'samtools view' format)
	elseif (ends_with($out_file, ".bam") && ends_with($reference_file, ".bam"))
	{
		$o = temp_file("_out.sam");
		$samtools_command = execApptainer("samtools", "samtools view", "$out_file", [$out_file], [], true);
		exec("{$samtools_command} | cut -f1-11 > $o 2>&1", $output, $return); //cut to ignore the read-group and other annotations
		$r = temp_file("_ref.sam");
		$samtools_command = execApptainer("samtools", "samtools view", "$reference_file", [$reference_file], [], true);
		exec("{$samtools_command} | cut -f1-11 > $r 2>&1", $output, $return); //cut to ignore the read-group and other annotations
		
		exec("diff -b $r $o > $logfile 2>&1", $output, $return);
		$passed = ($return==0);
	}
	//vcf
	elseif (ends_with($out_file, ".vcf") && ends_with($reference_file, ".vcf"))
	{
		$o = temp_file("_out.vcf");
		file_put_contents($o, implode("\n", load_vcf_normalized($out_file)));
		$r = temp_file("_ref.vcf");
		file_put_contents($r, implode("\n", load_vcf_normalized($reference_file)));
		
		exec("diff -b $r $o > $logfile 2>&1", $output, $return);
		$passed = ($return==0);
	}
	//diff
	else
	{
		$extras = "";
		if (!$compare_header_lines) $extras .= " -I '^[#@]'";
		exec("diff $extras $reference_file $out_file > $logfile 2>&1", $output, $return);

		$passed = ($return==0);
	}
	
	if ($passed)
	{
		$result = "PASSED";
		++$GLOBALS["passed"];
	}
	else
	{
		$result = "FAILED (see $logfile)";
		++$GLOBALS["failed"];
	}
	
	$bt = debug_backtrace();
	$caller = array_shift($bt);
	$file = basename($caller["file"]);
	$line = $caller["line"];
	print "  - $file:$line $result\n";
}

function check_file_tsv($out_file, $reference_file, $skip_comments = false, $skip_cols = "", $skip_comments_matching = "")
{
	check_test_started();
	
	$logfile = $out_file."_diff";
	
	$args = [];
	$args[] = "-in1 $reference_file";
	$args[] = "-in2 $out_file";
	if ($skip_comments) $args[] = "-skip_comments";
	if ($skip_cols!="") $args[] = "-skip_cols '$skip_cols'";
	if ($skip_comments_matching!="") $args[] = "-skip_comments_matching '$skip_comments_matching'";
	$tsvdiff_command = execApptainer("ngs-bits", "TsvDiff", implode(" ", $args), [$reference_file, $out_file], [], true);
	exec("$tsvdiff_command > $logfile 2>&1", $output, $return);
	$passed = ($return==0);
	
	if ($passed)
	{
		$result = "PASSED";
		++$GLOBALS["passed"];
	}
	else
	{
		$result = "FAILED (see $logfile)";
		++$GLOBALS["failed"];
	}
	
	$bt = debug_backtrace();
	$caller = array_shift($bt);
	$file = basename($caller["file"]);
	$line = $caller["line"];
	print "  - $file:$line $result\n";
}


/// Executes a command and checks that it does not return an error code
function check_exec($command, $fail = TRUE)
{
	$start_time = microtime(true);

	check_test_started();
	
	// execute command
	if ($GLOBALS["debug"]) print "  Executing: $command\n";
	exec($command." 2>&1", $output, $return);
	
	//check if output contains warning from PHP
	$warning = false;
	foreach ($output as $line)
	{
		if (starts_with($line, "PHP") && (contains($line, "Warning") || contains($line, "Error") || contains($line, "Notice")))
		{
			$warning = true;
		}
	}
	
	if (($return!=0 || $warning) && $fail)
	{
		++$GLOBALS["failed"];
		
		//write output to logfile
		$logfile = output_folder().random_string(4)."_output";
		file_put_contents($logfile, implode("\n", $output));
		
		//issue error message
		$bt = debug_backtrace();
		$caller = array_shift($bt);
		$file = basename($caller["file"]);
		$line = $caller["line"];
		print "  - $file:$line FAILED (see $logfile)\n";
	}
	else
	{
		$result = "PASSED";
		++$GLOBALS["passed"];

		$bt = debug_backtrace();
		$caller = array_shift($bt);
		$file = basename($caller["file"]);
		$line = $caller["line"];
		print "  - $file:$line $result\n";
	}
	
	if ($GLOBALS["debug"])
	{
		foreach($output as $line)
		{
			$line = trim($line);
			if ($line!="")
			{
				print "    output: ".$line."\n";
			}
		}
	}
	if ($GLOBALS["debug"])
	{
		print "    execution time: ".time_readable(microtime(true) - $start_time)."\n";
	}
	
	return $output;
}

/// Ends a test
function end_test()
{
	check_test_started();

	$failed = $GLOBALS["failed"]!=0;
	
	print ($failed ? "FAILED" : "FINISHED")." ".$GLOBALS["test_name"]." - passed: ".$GLOBALS["passed"]."/".($GLOBALS["passed"]+$GLOBALS["failed"])." - time: ".time_readable(microtime(true) - $GLOBALS["start_time"])."\n\n";
	
	if ($failed) exit(1);
}

///initializes NGSD if a SQL file with the correct name exists
function init_ngsd($tool_name)
{	
	//prepare arguments
	$args = [];
	$args[] = "-test";
	$in_files = [];
	
	//additional SQL commands to execute
	$sql_file = data_folder()."/{$tool_name}.sql";
	if (file_exists($sql_file))
	{
		$args[] = "-add {$sql_file}";
		$in_files[] = $sql_file;
	}
	
	//init
	list($stdout, $stderr, $exit_code) = execApptainer("ngs-bits", "NGSDInit", implode(" ", $args), $in_files);
}
?>
