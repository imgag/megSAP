<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function rtf_table_rows($rows, $col_widths, $text_aligns = array(), $border = false, $bold = false)
{
	$return = "";
	foreach($rows as $r)
	{
		if(count($col_widths)>0 && count($r)!=count($col_widths))	trigger_error("Less width information available than columns.",E_USER_ERROR);
		if(count($text_aligns)>0 && count($r)!=count($text_aligns))	trigger_error("Less align information available than columns.",E_USER_ERROR);

		$row = "{\trowd\trgaph70\fs20";		
		// column widths
		foreach($col_widths as $cw)
		{
			if($border)	$row .= "\clbrdrt\brdrw18\brdrs\clbrdrl\brdrw18\brdrs\clbrdrb\brdrw18\brdrs\clbrdrr\brdrw18\brdrs";
			$row .= "\cellx$cw";
		}
		// col data
		for($i=0;$i<count($r);++$i)
		{
			$c = $r[$i];
			$row .= "\pard\intbl\sa20\sb20";
			if(isset($text_aligns[$i]))	$row .= $text_aligns[$i]."";
			if($bold)	$row .= "\b";
			$row .= " {".$c."}\cell";
		}
		$row .= "\row}";
		
		$return .= "\n$row";
	}
	return $return;
}

function rtf_table_row($row, $col_widths, $text_aligns = array(), $border = false, $bold = false)
{
	return rtf_table_rows(array($row), $col_widths, $text_aligns, $border, $bold);
}


?>