<?php

require_once("functions.php");

/**
	@brief Base class for tools.
	
	This class support common tool tasks, i.e.:
	- Command line parsing
	- Logging
	- Execution of external tools
	- Error handling
	- Temporary file creation and cleanup
	
	@ingroup base
*/
class ToolBase
{
	private $name = "";
	private $description = "";
	private $params = array();
	
	private $log_file = null;
	private $log_id = "";
	private $start_time = null;
	private $debug = false;
	
	private $temp_files = array();
	private $temp_folders = array();
	private $error_occurred = false;
	
	/// Constructor
	function __construct($name, $description)
	{
        $this->name = $name;
        $this->description = $description;
		
		// check name
		$filename = basename($_SERVER['SCRIPT_FILENAME']);
		if ($filename != $name.".php")
		{			
			trigger_error("Tool name '$name' and file name '$filename' do not match!", E_USER_ERROR);
		}
		
		// Set start time
		$this->start_time = microtime(true);
		
		// Create random log_id. Overwritten if log_id command line flag is used
		$this->log_id = random_string(4);
		
		// register callbacks
		register_shutdown_function(array($this, 'shutdown'));
		set_error_handler(array($this, 'error_handler'));
	}
	
	/// Shutdown function (destructor cannot be used, because it is not called after 'trigger_error' is used)
	function shutdown()
	{
		//remove all temp files
		if (!$this->error_occurred && !$this->debug)
		{
			foreach($this->temp_files as $file)
			{
				if (file_exists($file) && !unlink($file))
				{
					trigger_error("Temporary file '$file' could not be deleted.", E_USER_NOTICE);
				}
			}
			
			foreach($this->temp_folders as $folder)
			{
			    if(is_dir($folder))
				{
					$this->removeTempFolder($folder);
				}
			}
		}
		
		//final log message
		$this->log("END ".$this->name);
		$this->log("Execution time of '".$this->name."': ".time_readable(microtime(true) - $this->start_time));
	}
	
	/// Remove temporary folder
	private function removeTempFolder($folder)
	{
		//only delete files within tmp folder
		if(!is_in_temp_folder($folder)) trigger_error("Folder '$folder' not within tmp-folder '$tmp'.", E_USER_ERROR);
		if (substr($folder, strlen($folder) - 1, 1)!='/') $folder .= '/';
		
		//get all files / sub_folders
		$files = scandir($folder);

		foreach ($files as $file) 
		{
			//exclude '.' and '..'
			if($file == "." || $file == "..") continue;

			if (is_dir($folder.$file)) 
			{
				//delete recursively
				$this->removeTempFolder($folder.$file);
			} 
			else 
			{
				//delete file
				unlink($folder.$file);
			}
		}
		
		//delete folder
		rmdir($folder);
	}
	
	/// Handles errors
	function error_handler($level, $message, $file, $line, $context="")
	{
		$type = "";
		//log error
		if($level & E_USER_ERROR)
		{
			$this->error_occurred = true;
			$type = "ERROR";
		}
		//log warning
		if($level & E_USER_WARNING)
		{
			$type = "WARNING";
		}
		//log error
		if($level & E_USER_NOTICE)
		{
			$type = "NOTICE";
		}
		
		//log and print known types
		if ($type!="")
		{
			if ($this->debug)
			{
				$tb = traceback();
				$this->log("$type: '$message' in $file:$line.", array($tb));
				$this->toStdErr("$type: '$message' in $file:$line.\n".$tb);
			}
			else
			{
				$this->log("$type: '$message' in $file:$line.");
				$this->toStdErr("$type: '$message' in $file:$line.");
			}
		}
		
		return false; //use PHP's error handler afterwards
	}

	/// @name Parameter handling members.
	//@{
	
	/// Adds a flag. It not used, the value parameter is set to 'false' (do not used isset!).
	function addFlag($name, $desc)
	{
		$this->params["-".$name] = array("flag", $desc, null, null, null, null);
	}
	
	/// Adds an input file parameter
	function addInfile($name, $desc, $opt, $check_readable = true)
	{
		$this->add($name, $desc, $opt, null, "infile", $check_readable);
	}

	/// Adds an input file array parameter
	function addInfileArray($name, $desc, $opt, $check_readable = true)
	{
		$this->add($name, $desc, $opt, null, "infile_array", $check_readable);
	}
		
	/// Adds an output file parameter. The variable is unset if optional and not given!
	function addOutfile($name, $desc, $opt)
	{
		$this->add($name, $desc, $opt, null, "outfile");
	}

	/// Adds an output file array parameter.
	function addOutfileArray($name, $desc, $opt)
	{
		$this->add($name, $desc, $opt, null, "outfile_array");
	}
	
	/// Adds an enum parameter
	function addEnum($name, $desc, $opt, $valid, $default = null)
	{
		if (isset($default) && !in_array($default, $valid))
		{
			trigger_error("Default value '$default' of enum '$name' not contained in valid values ".implode(",", $valid).".", E_USER_ERROR);
		}
		
		$this->add($name, $desc, $opt, $default, "enum", $valid);
	}

	/// Adds an string parameter
	function addString($name, $desc, $opt, $default = null)
	{
		$this->add($name, $desc, $opt, $default, "string");
	}

	/// Adds an string parameter
	function addStringArray($name, $desc, $opt, $default = null)
	{
		$this->add($name, $desc, $opt, $default, "string_array");
	}
	
	/// Adds an int parameter
	function addInt($name, $desc, $opt, $default = null)
	{
		$this->add($name, $desc, $opt, $default, "int");
	}

	/// Adds a float parameter
	function addFloat($name, $desc, $opt, $default = null)
	{
		$this->add($name, $desc, $opt, $default, "float");
	}

	/// Adds a parameter (internal)
	private function add($name, $desc, $opt, $default, $type, $option1 = null, $option2 = null)
	{
		//check for duplicate name
		if(isset($this->params["-".$name]))
		{
			trigger_error("Parameter name '-$name' given twice.", E_USER_ERROR);
		}
		
		if ($opt)
		{
			$this->params["-".$name] = array("opt", $desc, $type, $option1, $option2, $default);
		}
		else
		{
			if (isset($default))
			{
				trigger_error("Default value '$default' for non-optional parameter '$name' given.", E_USER_ERROR);
			}
			
			$this->params["-".$name] = array("par", $desc, $type, $option1, $option2, null);
		}
	}
	
	/// Main parsing function
	/// @note Use the PHP function extract() to import the result to the calling namespace
	function parse($argv)
	{	
		// Special paramters
		if (in_array("--help", $argv))
		{
			$this->printUsage();
			exit(0);
		}
		if(in_array("--version", $argv))
		{
			print "$this->name version ".$this->version()."\n";
			exit(0);
		}
		if(in_array("--debug", $argv))
		{
			$this->debug = true;
			$arg_index = array_search("--debug", $argv);
			array_splice($argv, $arg_index, 1, array());	
		}
		if(in_array("--log", $argv))
		{
			$arg_index = array_search("--log", $argv);
			$this->log_file = $argv[$arg_index+1];
			array_splice($argv, $arg_index, 2, array());	
		}
		if(in_array("--log_id", $argv))
		{
			$arg_index = array_search("--log_id", $argv);
			$this->log_id = $argv[$arg_index+1];
			array_splice($argv, $arg_index, 2, array());	
		}
		if(in_array("--conf", $argv))
		{
			$arg_index = array_search("--conf", $argv);
			$GLOBALS["path_ini"] = $argv[$arg_index+1];
			array_splice($argv, $arg_index, 2, array());	
		}
		
		//log tool start
		$this->log("START ".$this->name." (version: ".$this->version().")");
		$this->log("Debug: ".($this->debug ? "true" : "false"));
		
		//set default values
		$output = array();
		foreach($this->params as $name => $details)
		{
			if ($details[0] == "flag")
			{
				$output[substr($name,1)] = false;
			}
			else if ($details[0] == "opt" && (is_null($details[5]) || isset($details[5])))
			{
				$output[substr($name,1)] = $details[5];
			}
		}
		
		//parse command line
		$para_present = array();
		for ($i=1; $i<count($argv); ++$i)
		{
			$par = $argv[$i];
			if ($par[0] != "-")
			{
				$this->printUsage();
				trigger_error("Trailing argument '$par' given.", E_USER_ERROR);
			}
			
			
			if (!isset($this->params[$par]))
			{				
				$this->printUsage();
				trigger_error("Unknown parameter '$par' given.", E_USER_ERROR);
			}
			
			// handle parameters with arguments
			$para_type = $this->params[$par][0];
			if ($para_type=="par" || $para_type=="opt")
			{
				// get argument
				if (!isset($argv[$i + 1]) || ($argv[$i + 1][0]=="-" && isset($this->params[$argv[$i + 1]])))
				{
					$this->printUsage();
					trigger_error("Parameter '$par' given without argument.", E_USER_ERROR);
				}
				$arg = $argv[$i + 1];
				
				if (in_array($par, $para_present))
				{
					$this->printUsage();
					trigger_error("Parameter '$par' given more than once.", E_USER_ERROR);
				}
				$para_present[] = $par;
				
				// check argument data types
				$argu_type = $this->params[$par][2];
				$argu_rest = $this->params[$par][3];
				if ($argu_type == "infile")
				{
					$check_readable = $this->params[$par][3];
					if (realpath($arg) === false)
					{
						$this->printUsage();
						trigger_error("Could not find realpath for input file '".$argv[$i+1]."' given for list parameter '$par'. Keeping relative path.", E_USER_WARNING);
					}
					else
					{
						$arg = realpath($arg);
					}
					
					if ($check_readable && !is_readable($arg))
					{
						$this->printUsage();
						trigger_error("Input file '$arg' given for parameter '$par' is not readable.", E_USER_ERROR);
					}
				}
				else if($argu_type == "infile_array")
				{
					$arg = array();
					while(isset($argv[$i + 1]) && $argv[$i + 1][0]!="-")
					{
						//print $argv[$i+1]."\n";
						if (realpath($argv[$i+1]) === false)
						{
							trigger_error("Could not find realpath for input file '".$argv[$i+1]."' given for list parameter '$par'. Keeping relative path.", E_USER_WARNING);
							$arg[] = $argv[$i+1];
						}
						else
						{
							$arg[] = realpath($argv[$i+1]);
						}
						
						$check_readable = $this->params[$par][3];
						if ($check_readable && !is_readable(end($arg)))
						{
							$this->printUsage();
							trigger_error("Input file '".end($arg)."' given for list parameter '$par' is not readable.", E_USER_ERROR);
						}
						++$i;
					}
					--$i;
				}
				else if ($argu_type == "outfile")
				{
					if (!is_writable2($arg))
					{
						$this->printUsage();
						trigger_error("Output file '$arg' given for parameter '$par' is not writable.", E_USER_ERROR);
					}
				}
				else if($argu_type == "outfile_array")
				{
					$arg = array();
					while(isset($argv[$i + 1]) && $argv[$i + 1][0]!="-")
					{
						$arg[] = $argv[$i+1];
						if (!is_writable2(end($arg)))
						{
							$this->printUsage();
							trigger_error("Output file '".$argv[$i+1]."' given for list parameter '$par' is not writable.", E_USER_ERROR);
						}
						++$i;
					}
					--$i;
				}
				else if ($argu_type == "enum")
				{
					if (!in_array($arg, $argu_rest))
					{
						$this->printUsage();
						trigger_error("Invalid value '$arg' given for restriced parameter '$par'.", E_USER_ERROR);
					}
				}
				else if ($argu_type == "string")
				{
					// nothing to check
				}
				else if ($argu_type == "string_array")
				{
					$arg = array();
					
					while(isset($argv[$i + 1]) && $argv[$i + 1][0]!="-")
					{
						$arg[] = $argv[$i+1];
						++$i;
					}
					--$i;
				}
				else if ($argu_type == "int")
				{
					if ($arg != (string)(int)$arg)
					{
						$this->printUsage();
						trigger_error("Invalid value '$arg' given for integer parameter '$par'.", E_USER_ERROR);
					}
				}
				else if ($argu_type == "float")
				{
					if ($arg != (string)(float)$arg)
					{
						$this->printUsage();
						trigger_error("Invalid value '$arg' given for float parameter '$par'.", E_USER_ERROR);
					}
				}
				
				$output[substr($par,1)] = $arg;
				
				++$i;
			}
			// handle flags
			else if ($para_type=="flag")
			{
				$output[substr($par,1)] = true;
			}
		}
		
		// check missing parameters
		foreach($this->params as $name => $details)
		{
			if ($details[0] == "par" && !in_array($name, $para_present))
			{
				$this->printUsage();
				trigger_error("Mandatory parameter '$name' not given.", E_USER_ERROR);
			}
		}
		
		//log parameters
		if (count($output)!=0)
		{
			$length = max(array_map('strlen', array_keys($output)));
			foreach($output as $name => $value)
			{
				if(is_array($value))
				{
					$value = implode(" ", $value);
				}
				$this->log("Parameter: ".str_pad($name, $length, " ")." = $value");
			}
		}
		
		return $output;
	}
	
	/// Prints the usage information
	function printUsage()
	{
		// find out longest parameter and agrument name
		$max_name = 0;
		$max_arg = 0;
		foreach($this->params as $name => $details)
		{
			$max_name = max($max_name, strlen($name));
			$max_arg = max($max_arg, strlen($details[2]));
		}
		$offset = max(18, $max_name + $max_arg + 6);
		$indent = "  ".str_pad("", $offset, " ");
		
		print $this->name." (version ".$this->version().")\n";
		print "\n";
		print $this->description."\n";
		
		// Mandatory parameters
		$mandatory = array();
		foreach($this->params as $name => $details)
		{
			if ($details[0]=="par")
			{
				$mandatory[] = "  ".str_pad($name." <".$details[2].">", $offset, " ").str_replace("\n", "\n$indent", $details[1])."\n";
				
				if ($details[2]=="enum")
				{
					$mandatory[] = $indent."Valid are: '".implode("', '", $details[3])."'.\n";
				}
			}
		}
		if (count($mandatory))
		{
			print "\n";
			print "Mandatory parameters:\n";
			print implode("", $mandatory);
		}

		// Optional parameters
		$optional = array();
		foreach($this->params as $name => $details)
		{
			if ($details[0]=="opt")
			{
				$optional[] = "  ".str_pad($name." <".$details[2].">", $offset, " ").str_replace("\n", "\n$indent", $details[1])."\n";
				
				if ($details[2]=="enum")
				{
					$optional[] = $indent."Valid are: '".implode("', '", $details[3])."'.\n";
				}
				
				if (isset($details[5]))
				{
					$optional[] = $indent."Default is: '".$details[5]."'.\n";
				}
				
			}
			else if ($details[0] == "flag")
			{
				$optional[] = "  ".str_pad($name, $offset, " ").str_replace("\n", "\n$indent", $details[1])."\n";
			}
		}
		if (count($optional))
		{
			print "\n";
			print "Optional parameters:\n";
			print implode("", $optional);
		}

		print "\n";
		print "Special parameters:\n";
		print "  ".str_pad("--help", $offset, " ")."Shows this help.\n";
		print "  ".str_pad("--version", $offset, " ")."Prints version and exits.\n";
		print "  ".str_pad("--debug", $offset, " ")."Enables verbose error messages including traceback and disables automatic deletion of temporary files/folders.\n";
		print "  ".str_pad("--log <file>", $offset, " ")."Logs to the specified file. Use '-' to log to STDOUT.\n";
		print "  ".str_pad("--log_id <id>", $offset, " ")."Uses the given identifier for logging.\n";
		print "  ".str_pad("--conf <file>", $offset, " ")."Uses the given INI file instead of the default settings.ini.\n";
		print "  ".str_pad("--tdx", $offset, " ")."Writes a Tool Defition XML file.\n";
		
		print "\n";
		print "\n";
	}

	//@}
	
	
	/// @name Logging members
	//@{
	
	/// Logs a message and optional additional info.
	public function log($message, $add_info = array())
	{
		if (!isset($this->log_file)) return;
		
		//skip excessive output of some tools
		$skip_substrings = [
							"[M::mem_pestat]", //BWA mem
							"[M::process]", //BWA mem
							"[M::mem_process_seqs]", //BWA mem
							"PROCESS_REGION_MSECS", //ABRA2
							"Processing chromosome chunk:", //ABRA2
							"Clock time in Chromosome:", //ABRA2
							"WARNING(freebayes): Could not find any mapped reads in target region", //freebayes
							"Subroutine Bio::DB::IndexedBase::_strip_crnl redefined", //VEP
							"No training configuration", //SpliceAI
							"ms/step", //SpliceAI
						   ];
		
		if ($this->log_file=="-") // to STDOUT
		{
			print $message."\n";
			foreach($add_info as $line)
			{
				if (trim($line)=="") continue;
				
				foreach($skip_substrings as $substring)
				{
					if (contains($line, $substring)) continue(2);
				}
				
				print "  ".trim($line)."\n";
			}
		}
		else // to file
		{
			$prefix = get_timestamp()."\t".$this->log_id."\t";
			
			$lines = array();
			$lines[] = $prefix."$message\n";
			foreach($add_info as $line)
			{
				if (trim($line)=="") continue;

				foreach($skip_substrings as $substring)
				{
					if (contains($line, $substring)) continue(2);
				}
				
				$lines[] = $prefix."    ".strtr($line, array("\n" => "", "\c" => "", "\t" => "  "))."\n";
			}
			
			$result = file_put_contents($this->log_file, $lines, LOCK_EX | FILE_APPEND);
			if ($result===FALSE)
			{
				trigger_error("Could not write to ".$this->log_file."!", E_USER_ERROR);
			}
		}
	}
	
	
	///Logs server name, user name and working directory.
	public function logServerEnvronment()
	{
		$server = trim(implode(" ", exec2("hostname -f")[0]));
		$this->log("Environment/server: {$server}");
		$os = trim(implode(" ", exec2("cat /etc/issue")[0]));
		$this->log("Environment/OS: {$os}");
		$user = trim(exec('whoami'));
		$this->log("Environment/user: {$user}");
		$lang = trim(exec('echo $LANG'));
		$this->log("Environment/LANG: {$lang}");
		$this->log("Environment/container platform: ".container_platform(true));
	}
	
	/// Prints a message (string or array of strings) to the stderr stream
	function toStderr($message)
	{
		if (!is_array($message))
		{
			$message = array($message);
		}
		
		$handle = fopen2('php://stderr', 'w');
		foreach($message as $line)
		{
			fwrite($handle,nl_trim($line)."\n");
		}
		fclose($handle);
	}
	
	/// Prints/logs a warning message.
	public function warning($message, $info = array())
	{
		$message = "WARNING: ".$message;

		//print infos
		print $message."\n";
		
		//log error
		$this->log($message, $info);
	}

	//@}
	
	//tries to extract the version number from a command
	function extractVersion($command)
	{
		$version = "n/a";
		
		//tools with --version argument
		$output = array();
		exec($command." --version 2>&1", $output, $return);
		if ($return==0 && preg_match("/[0-9]+[\.-_][0-9\.]+[0-9A-Za-z_-]*$/", $output[0], $hits))
		{
			$version = $hits[0];
		}
		
		//tools with -version argument
		if ($version=="n/a")
		{
			$output = array();
			exec($command." -version 2>&1", $output, $return);
			if ($return==0 && preg_match("/[0-9]+[\.-][0-9\.]+[0-9A-Za-z_-]*$/", $output[0], $hits))
			{
				$version = $hits[0];
			}
		}
		
		return $version;
	}

	/**
		@brief Executes the command and returns an array with STDOUT, STDERR and exit code.
	*/
	function exec($command, $parameters, $log_output=true, $abort_on_error=true, $warn_on_error=true)
	{
		if (is_array($command) || is_array($parameters))
		{
			print_r($command);
			print_r($parameters);
			die;
		}
		//prevent execution of pipes - exit code is not handled correctly with pipes!
		$command_and_parameters = $command." ".$parameters;
		if(contains($command_and_parameters, "|"))
		{
			trigger_error("Error in 'exec' method call: Command must not contain pipe symbol '|'. Please use the execPipeline function instead!\n$command_and_parameters", E_USER_ERROR);
		}
		
		//log call
		if($log_output)
		{
			$add_info = array();
			$add_info[] = "version    = ".$this->extractVersion($command);
			$add_info[] = "parameters = $parameters";
			$this->log("Calling external tool '$command'", $add_info);
		}
		
		$pid = getmypid();
		//execute call - pipe stdout/stderr to file
		$stdout_file = $this->tempFile(".stdout", "megSAP_exec_pid{$pid}_");
		$stderr_file = $this->tempFile(".stderr", "megSAP_exec_pid{$pid}_");
		$exec_start = microtime(true);
		$proc = proc_open($command_and_parameters, array(1 => array('file',$stdout_file,'w'), 2 => array('file',$stderr_file,'w')), $pipes);
		if ($proc===false) trigger_error("Could not start process with ToolBase::exec function!\nCommand: {$command}\nParameters: {$parameters}", E_USER_ERROR);
		
		//get stdout, stderr and exit code
		$return = proc_close($proc);
		$stdout = explode("\n", rtrim(file_get_contents($stdout_file)));
		$stderr = explode("\n", rtrim(file_get_contents($stderr_file)));
		$this->deleteTempFile($stdout_file);
		$this->deleteTempFile($stderr_file);
			
		//log relevant information
		if (($log_output || $return!=0) && count($stdout)>0)
		{
			$this->log("Stdout of '$command':", $stdout);
		}
		if (($log_output || $return!=0) && count($stderr)>0)
		{
			$this->log("Stderr of '$command':", $stderr);
		}
		if ($log_output)
		{
			$this->log("Execution time of '$command': ".time_readable(microtime(true) - $exec_start));
		}
		
		//abort/warn if failed
		if ($return!=0 && ($abort_on_error || $warn_on_error))
		{
			$this->toStderr($stdout);
			$this->toStderr($stderr);
			trigger_error("Call of external tool '$command' returned error code '$return'.", $abort_on_error ? E_USER_ERROR : E_USER_WARNING);
		}
		
		return array($stdout, $stderr, $return);
	}

	/**
		@brief Executes a pipeline (array of commands and parameters) and returns an array with STDOUT and STDERR.
		
		If the call exits with an error code, further execution of the calling script is aborted.
		If no version number is given, an automatic extraction of the version number is tried.
	*/
	function execPipeline($commands_params, $name, $log_output=true)
	{
		$add_info = array();
		$parts = array();
		$stderr_files = [];
		foreach($commands_params as $call)
		{
			list($command, $params) = $call;
			
			$stderr_file = $this->tempFile(".stderr");
			$parts[] = "{$command} {$params} 2>{$stderr_file}";
			$stderr_files[] = $stderr_file;
			$add_info[] = "command ".count($parts)."  = $command";
			if (!isset($call[2]) || $call[2]==true)
			{
				$add_info[] = "version    = ".$this->extractVersion($command);
			}
			$add_info[] = "parameters = $params";	
		}
			
		//log call
		$pipeline = implode(" | ", $parts);		
		$add_info[] = "pipeline   = $pipeline";			
		if($log_output)
		{
			$this->log("Calling '$name' pipeline", $add_info);
		}
			
		$exec_start = microtime(true);
		exec("bash -c \"set -o pipefail && ".strtr($pipeline, array("\""=>"\\\""))."\"", $stdout, $return); //pipefail is needed because otherwise only the exit code of the last tool is returned!
		
		//log stdout
		if (($log_output || $return!=0) && count($stdout)>0)
		{
			$this->log("Stdout of '$name' pipeline:", $stdout);
		}
		//log stderr
		$stderr = [];
		if($log_output || $return!=0)
		{
			foreach($stderr_files as $stderr_file)
			{
				foreach(file($stderr_file) as $line)
				{
					if (stripos($line, "WARNING: Error changing the container working directory") !== false) continue;
					if (trim($line)!="") $stderr[] = $line;
				}
			}
			if (count($stderr)>0)
			{
				$this->log("Stderr of '$name' pipeline:", $stderr);
			}
		}
		//log execution time
		if ($log_output)
		{
			$this->log("Execution time of '$name' pipeline: ".time_readable(microtime(true) - $exec_start));
		}
		
		//abort on error
		foreach($stderr as $line)
		{
			if ($return!=0) break;
		}
		if ($return!=0)
		{	
			$this->toStderr($stderr);
			trigger_error("Call of '$name' pipeline returned error code '$return'.", E_USER_ERROR);
		}
		
		//return results
		return array($stdout, $stderr);
	}

	/**
	 	@brief Executes a command inside a given Apptainer container and returns an array with STDOUT, STDERR and exit code.
	 */
	function execApptainer($container, $command, $parameters, $in_files=[], $out_files=[], $command_only=false, $log_output=true, $abort_on_error=true, $warn_on_error=true, $gpu_container=false)
	{
		//execute local ngs-bits installation if ngs-bits_local path is set in settings.ini
		$ngsbits_local = get_path("ngs-bits_local", false);
		if ($container == "ngs-bits" && $ngsbits_local != "")
		{
			if ($command_only)
			{
				$ngsbits_command = execApptainer($container, $command, $parameters, $in_files, $out_files, true, false, $abort_on_error);
				return $ngsbits_command;
			}
			else
			{
				list($stdout, $stderr, $return) = execApptainer($container, $command, $parameters, $in_files, $out_files, false, false, $abort_on_error);
				return array($stdout, $stderr, $return);
			}
		}

		//get apptainer command, bind paths and container path from execApptainer function in functions.php
		list($apptainer_command, $bind_paths, $container_path, $container_version) = execApptainer($container, $command, $parameters, $in_files, $out_files, false, true, true, $gpu_container);
		
		//if command only option is true, only the apptainer command is being return, without execution
		if($command_only) 
		{
			//$this->log("DEBUG: Apptainer command:\t", array($apptainer_command));
			return $apptainer_command;
		}
		//log call
		if($log_output)
		{
			$version_command = "singularity exec $container_path $command";
			$add_info = array();
			foreach($in_files as $in_file)
			{
				$add_info[] = "input file          = ".$in_file;
			}
			foreach($out_files as $out_file)
			{
				$add_info[] = "output file         = ".$out_file;
			}
			foreach($bind_paths as $bind_path)
			{
				$add_info[] = "bind path           = ".$bind_path;
			}
			$add_info[] = "container path      = ".$container_path;
			$add_info[] = "tool version        = ".$this->extractVersion($version_command);
			$add_info[] = "parameters          = $parameters";
			$this->log("Calling external tool '$command' in container '".basename2($container_path)."'", $add_info);
		}
				
/* 		$this->log("DEBUG: Apptainer version: ".$this->extractVersion("apptainer"));
		$this->log("DEBUG: Apptainer command: ".$apptainer_command);  */
		
		
		$pid = getmypid();
		//execute call - pipe stdout/stderr to file
		$stdout_file = $this->tempFile(".stdout", "megSAP_exec_pid{$pid}_");
		$stderr_file = $this->tempFile(".stderr", "megSAP_exec_pid{$pid}_");
		$exec_start = microtime(true);
		$proc = proc_open($apptainer_command, array(1 => array('file',$stdout_file,'w'), 2 => array('file',$stderr_file,'w')), $pipes);
		if ($proc===false) trigger_error("Could not start process with ToolBase::execApptainer function!\nContainer: {$container_path}\nCommand: {$command}\nParameters: {$parameters}", E_USER_ERROR);
		
		//get stdout, stderr and exit code
		$return = proc_close($proc);
		$stdout = explode("\n", rtrim(file_get_contents($stdout_file)));
		$stderr = explode("\n", rtrim(file_get_contents($stderr_file)));
		$this->deleteTempFile($stdout_file);
		$this->deleteTempFile($stderr_file);
			
		// Clean stderr by removing specific warning lines
		$filtered_stderr_lines = array_filter($stderr, function($line) 
		{
			return stripos($line, "WARNING: Error changing the container working directory") === false;
		});

		//log relevant information
		if (($log_output || $return!=0) && count($stdout)>0)
		{
			$this->log("Stdout of '$command':", $stdout);
		}
		if (($log_output || $return!=0) && count($filtered_stderr_lines)>0)
		{
			$this->log("Stderr of '$command':", $filtered_stderr_lines);
		}
		if ($log_output)
		{
			$this->log("Execution time of '$command': ".time_readable(microtime(true) - $exec_start));
		}
		
		//abort/warn if failed
		if ($return!=0 && ($abort_on_error || $warn_on_error))
		{
			$this->toStderr($stdout);
			$this->toStderr($filtered_stderr_lines);
			trigger_error("Call of external tool '$command' returned error code '$return'.", $abort_on_error ? E_USER_ERROR : E_USER_WARNING);
		}
		
		return array($stdout, $filtered_stderr_lines, $return);
	}
	//TODO Marc/Leon: the exit code seems to be 0 in case of error. See /mnt/storage2/users/ahsturm1/scripts/2025_06_18_missing_variants_in_trio/ and TARGETED workaround in merge_gvcf.php
	/**
	 	@brief Executes a list of commands in parallel
	*/
	function execParallel($command_list, $threads, $log_output=true, $abort_on_error=true, $log_stdout=true, $log_summary=false)
	{
		
		$running = array();

		//create array with return values
		$return  = array();
		$return["jobs"] = array();
		$return["jobs_failed"] = array();
		$return["jobs_successfully_finished"] = array();
		$return["n_jobs_started"] = 0;
		$return["n_jobs_successfully_finished"] = 0;
		$return["n_jobs_failed"] = 0;

		while (true) 
		{
			//wait 50ms
			usleep(50000);
			
			//jobs have been processed > exit
			if (count($command_list)==0 && count($running)==0) break;
			
			//start new sub-processed while there are threads unused
			while(count($running)<$threads && count($command_list)>0)
			{
				list($name, $command) = array_shift($command_list);
				list($command, $parameters) = explode(" ", $command, 2);
				$job_id = $name."_".random_string(8);
				$return["jobs"][$job_id] = array(); 

				$descriptors = [
					1 => ['pipe', 'w'], // stdout
					2 => ['pipe', 'w'], // stderr
				];
				$process = proc_open($command.' '.$parameters, $descriptors, $pipes);
				if (is_resource($process)) 
				{
					$status = proc_get_status($process);
					//log start
					$add_info = array();
					$add_info[] = "version    = ".$this->extractVersion($command);
					$add_info[] = "parameters = {$parameters}";
					$add_info[] = "pid = ".$status['pid'];
					if($log_output) $this->log("Executing command in background '{$command}'", $add_info);
					$running[$job_id] = [
						"process" => $process,
						"pipes" => $pipes,
						"pid" => $status['pid'],
						"cmd" => $command.' '.$parameters,
						"start_time" => microtime(true),
					];

					//log in return array
					$return["jobs"][$job_id]["command"] = $command;
					$return["jobs"][$job_id]["parameters"] = $parameters;
					$return["jobs"][$job_id]["pid"] = $status['pid'];
					$return["n_jobs_started"]++;
				}
				else
				{
					trigger_error("Couldn't start process for command '{$command}'!", E_USER_ERROR);
				}
			}

			//check if running jobs have finished
			foreach ($running as $job_id => $info) 
			{
				$status = proc_get_status($info["process"]);
				if (!$status["running"])
				{
					$job_aborted = false;
					
					// get info about this finished job
					$stdout = stream_get_contents($info['pipes'][1]);
    				$stderr = stream_get_contents($info['pipes'][2]);
					fclose($info['pipes'][1]);
					fclose($info['pipes'][2]);
					$exit_code = proc_close($info["process"]);

					//logging
					$add_info = array();
					if($log_stdout)
					{
						if ($stdout!="")
						{
							$add_info[] = "STDOUT:";
							foreach(explode("\n", $stdout) as $line)
							{
								$add_info[] = nl_trim($line);
							}
						}
					}
					if ($stderr!="")
					{
						$add_info[] = "STDERR:";
						foreach(explode("\n", $stderr) as $line)
						{
							$add_info[] = nl_trim($line);
						}
					}
					$add_info[] = "EXIT CODE: ".$exit_code;
					if(!is_numeric($exit_code) || $exit_code!=0)
					{
						$job_aborted = true;
					}

					//log in return array
					$return["jobs"][$job_id]["stdout"] = explode("\n", $stdout);
					$return["jobs"][$job_id]["stderr"] = explode("\n", $stderr);
					$return["jobs"][$job_id]["exit_code"] = $exit_code;
					if ($job_aborted)
					{
						$return["jobs_failed"][] = $job_id;
						$return["n_jobs_failed"]++;
					}
					else
					{
						$return["n_jobs_successfully_finished"]++;
					}
					
					
					//log output
					if($log_output) $this->log("Finshed processing job {$job_id} in ".time_readable(microtime(true)-$info["start_time"]), $add_info);
					
					//abort if failed
					if ($job_aborted) trigger_error("Processing of job {$job_id} failed: ".$stderr, $abort_on_error ? E_USER_ERROR : E_USER_WARNING);
					
					unset($running[$job_id]);

				}
				//else: job is still running -> continue
			}	
		}

		if ($log_summary)
		{
			foreach ($return["jobs"] as $job_id => $info) 
			{
				if ($info["exit_code"] == 0)
				{
					$this->log("Job {$job_id} finished successfully (exit code: 0)");
					$this->log("\t command:\t\t".$info["command"] );
					$this->log("\t parameters:\t".$info["parameters"] );
				}
				else
				{
					$this->log("Job {$job_id} failed with exit code ".$info["exit_code"]);
					$this->log("\t command:\t\t".$info["command"] );
					$this->log("\t parameters:\t".$info["parameters"] );
					$this->log("\t stdout:\t\t".implode("\n\t\t\t\t\t", $info["stdout"]));
					$this->log("\t stderr:\t\t".implode("\n\t\t\t\t\t", $info["stderr"]));
				}
			}

			$this->log("SUMMARY: {$return['n_jobs_successfully_finished']} of {$return['n_jobs_started']} jobs finished successfully. ({$return['n_jobs_failed']} failed)");

		return $return;
		}
	}
	
	/// Returns the version of this tool
	function version()
	{
		return repository_revision();
	}
	
	/**
		@brief Executes a other tool which is based on ToolBase and returns an array with the output to STDOUT and STDERR.
		
		The output of the tool to stdout and stderr is logged. The version number is automatically determined and logged.
		
		If the call exits with an error code, further execution of the calling script is aborted.
	*/
	function execTool($command, $parameters, $abort_on_error=true)
	{
		//prepend php and path
		$command = "php ".repository_basedir()."/src/".$command;
		
		//determine tool version
		$output = array();
		exec($command." --version", $output, $return);
		preg_match("/[A-Za-z0-9\.-]+$/", $output[0], $hits);
		$version = $hits[0];
		
		//create call id
		$call_id = random_string(4);

		//log call
		$add_info = array();
		$add_info[] = "version    = $version";
		$add_info[] = "parameters = $parameters";
		$add_info[] = "log_id     = $call_id";
		$this->log("Calling external tool '$command'", $add_info);
		
		//execute call and pipe stderr stream to file
		$temp = $this->tempFile(".stderr");
		$addinfo = "";
		if (isset($this->log_file) && !contains($parameters, "--log"))
		{
			$addinfo .= "--log $this->log_file --log_id $call_id ";
		}
		if (isset($GLOBALS["path_ini"]) && !contains($parameters, "--conf"))
		{
			$addinfo .= "--conf ".$GLOBALS["path_ini"]." ";
		}
		if ($this->debug && !contains($parameters, "--debug"))
		{
			$addinfo .= "--debug ";
		}
		
		exec($command." ".$parameters." $addinfo 2>$temp", $stdout, $return);
		//print($command." ".$parameters." $addinfo 2>$temp");
		$stderr = file($temp);
		
		//log output streams
		if (count($stdout)>0)
		{
			$this->log("Stdout of '$command':", $stdout);
		}
		if (count($stderr)>0)
		{
			$this->log("Stderr of '$command':", $stderr);
		}
		
		//abort on error
		if ($abort_on_error && $return != 0)
		{	
			$this->toStderr($stderr);
			trigger_error("Call of external tool '$command' returned error code '$return'.", $abort_on_error ? E_USER_ERROR : E_USER_WARNING);
		}
		
		//return results
		return array($stdout, $stderr, $return);
	}
	
	/**
		@brief Returns a temporary file name
		
		@param suffix File suffix to append.
		@param prefix File prefix to prepend (default is the tool name).
		
		@note The file is deleted in the destructor of this class.
	*/
	function tempFile($suffix = "", $prefix = null, $directory = null)
	{
		if (!isset($prefix))
		{
			$prefix = $this->name."_";
		}
		$file = temp_file($suffix, $prefix, $directory);
		$this->temp_files[] = $file;
		
		return $file;
	}
	
	/**
		@brief Returns a temporary folder name
		
		@param prefix File prefix to prepend (default is the tool name).
		
		@note The folder is deleted in the destructor of this class.
	*/
	function tempFolder($prefix = null)
	{
		if (!isset($prefix))
		{
			$prefix = $this->name."_";
		}
		$folder = temp_folder($prefix);
		$this->temp_folders[] = $folder;
		
		return $folder;
	}	
	
	/**
		@brief Adds a file to temporary-files.
		
		@note The file is deleted in the destructor of this class.
	*/
	function addTempFile($file)
	{	
		if(!file_exists($file))
		{
			trigger_error("Temporary file '$file' does not exist!", E_USER_ERROR);
		}

		$this->temp_files[] = $file;
	}
	

	/**
		@brief Immediately deletes a temporary file, instead of deleting it in the destructor
	*/
	function deleteTempFile($file)
	{	
		if(!in_array($file, $this->temp_files))
		{
			trigger_error("Temporary file '$file' for which deletion was requested is unknown!", E_USER_ERROR);
		}
		
		if (file_exists($file) && !unlink($file))
		{
			trigger_error("Temporary file '$file' could not be deleted.", E_USER_NOTICE);
		}
		
		//remove temp file from list
		$new_temp_files = array();
		foreach($this->temp_files as $temp)
		{
			if ($temp==$file) continue;
			$new_temp_files[] = $temp;
		}
		$this->temp_files = $new_temp_files;
	}
	
	/**
		@brief Returns the log-file
		
	*/
	function getLogFile()
	{	
		return $this->log_file;
	}
	
	/**
		@brief Sets the log-file
		
	*/
	function setLogFile($file)
	{	
		$this->log_file = $file;
	}
	
	///Sort BAM file (optional build parameter to speed-up cram sorting)
	function sortBam($in, $out, $threads, $build, $by_name = FALSE)
	{	
		$threads -= 1; //number of additional threads, that's why -1
		$tmp_for_sorting = $this->tempFile();
		$this->execApptainer("samtools", "samtools sort", "-T {$tmp_for_sorting} -@ {$threads}".($by_name?" -n":"")." -m 1G --reference ".genome_fasta($build)." -o $out $in", [$in, genome_fasta($build)], [dirname($out)]);
	}
	
	///Index BAM file
	function indexBam($bam, $threads)
	{
		$this->execApptainer("samtools", "samtools index", "-@ {$threads} $bam", [$bam]);
	}
	
	///Move file with error checks
	function moveFile($from, $to)
	{
		$start = microtime(true);
		
		if (file_exists($to) && !unlink($to))
		{
			trigger_error("Could not move '$from' to '$to': Could not remove target!", E_USER_ERROR);
		}
		
		$this->exec("mv", "$from $to", false);
		
		$this->log("Execution time of moving '$from' to '$to': ".time_readable(microtime(true) - $start));		
	}
	
	///Copy file with error checks
	function copyFile($from, $to)
	{
		$start = microtime(true);
		
		if (file_exists($to) && !unlink($to))
		{
			trigger_error("Could not copy '$from' to '$to': Could not remove target!", E_USER_ERROR);
		}
		
		$this->exec("cp", "$from $to", false);
		
		$this->log("Execution time of copying '$from' to '$to': ".time_readable(microtime(true) - $start));		
	}
	
}

?>