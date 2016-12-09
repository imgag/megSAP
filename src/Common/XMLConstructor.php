<?php

require_once("functions.php");

/**
	@brief Convenient construction of an XML file (including validation).
	
	@ingroup xml
*/
class XMLConstructor
{
	private $output = array();
	private $xml_tags = array();
	
	///Constructor
	function __construct()
	{
        $this->output[] = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
		$this->xml_tags = array();
    }
	
	/// Opens an XML tag
	function openTag($name, $attr = array(), $cdata = "", $close=false)
	{
		$text = str_repeat("\t", count($this->xml_tags))."<".$name;
		foreach($attr as $n => $v)
		{
			$text .= " ".$n."=\"".$v."\"";
		}
		
		if ($cdata=="" && $close)
		{
			$text .= "/>";
		}
		else if ($cdata!="" && $close)
		{
			$text .= ">".htmlspecialchars($cdata)."</$name>";
		}
		else
		{
			$text .= ">".htmlspecialchars($cdata);
		}
		
		$this->output[] = $text."\n";
		
		if (!$close)
		{
			$this->xml_tags[] = $name;
		}
	}

	///Closes an open XML tag
	function closeTag($check_name = "")
	{
		if (count($this->xml_tags)==0)
		{
			trigger_error("Internal error: Closing XML tag requested, but all tags are already closed.", E_USER_ERROR);
		}

		$name = array_pop($this->xml_tags);
		
		if ($check_name!="" && $check_name!=$name)
		{
			trigger_error("Internal error: Cannot close XML tag '$check_name' beause the last opened tag is '$name'.", E_USER_ERROR);
		}
		
		$this->output[] = str_repeat("\t", count($this->xml_tags))."</".$name.">\n";
	}
	
	///Adds (opens and closes) an XML tag
	function AddTag($name, $attr = array(), $cdata = "")
	{
		$this->OpenTag($name, $attr, $cdata, true);
	}

	///Adds CDATA to an XML file
	function addCDATA($cdata)
	{
		$index = count($this->output)-1;
		$this->output[$index] .= htmlspecialchars($cdata);
	}
	
	///Returns if the output is wellformed
	function isWellformed($print_errors = true)
	{
		$messages = array();
		if (!xml_is_wellformed(implode("", $this->Output()), $messages))
		{
			if($print_errors)
			{
				xml_print_messages($messages);
			}
			
			return false;
		}
		
		return true;
	}
	
	///Returns the output string
	function output()
	{
		return $this->output;
	}
}


?>