<?php

set_error_handler("basic_error_handler");
declare(ticks = 1);

//functions
require_once("functions.php");
require_once("genomics.php");

//classes
require_once("XMLConstructor.php");
require_once("ToolBase.php");
require_once("Matrix.php");
require_once("Obo.php");
require_once("db.php");
require_once("rtf.php");

//basic error_handler for all libraries, reports only E_USER_ERROR
function basic_error_handler($level, $message, $file, $line, $context)
{
        if ($level & E_USER_ERROR)
        {
                if (!is_array($message))
		{
			$message = array($message);
		}
		
		$handle = fopen('php://stderr', 'w');
		foreach($message as $line)
		{
			fwrite($handle,"ERROR: ".nl_trim($line)."\n");
		}
		fclose($handle);
        }
        return false; //use PHP's error handler afterwards
}
?>