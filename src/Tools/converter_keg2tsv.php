<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

#parse command line arguments
$parser = new ToolBase("converter_keg2tsv", "\$Rev: 309 $", "");
$parser->addInfile("in", "Convert keg-format (exports from KEGG database) to tsv-Format.", false);
$parser->addOutfile("out", "TSV-file containing gene names.", false);
extract($parser->parse($argv));

$in_file = file($in, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);

//extract data
$i=0;
$rows = array();
while($in_file[$i]!="!")	++$i;
while($in_file[++$i]!="!")
{
	$rows[] = $in_file[$i];
}

//sort data
$data = array();
$sections = array();
$has_elements = false;
for($i=0;$i<count($rows);++$i)
{
	if($rows[$i][0] == "#")
	{
		//loop through prev_section, nb.: move tree one level up (before A)
		$rows[$i] = chr($prev_section);
		continue;
	}
	
	$current_section = ord(substr($rows[$i],0,1));
	$row = strip_tags(trim(substr($rows[$i],1)));
	
	if($i>0)	$prev_section = ord($rows[$i-1]);
	else	$prev_section = ord("A");

	if($i<count($rows))	$next_section = ord("A");
	else	$next_section = ord($rows[$i+1]);

	//sections
	if(!is_numeric($row[0]))
	{
		$length = $current_section-ord("A");
		
		if($prev_section < $current_section)
		{
			//simply add element
			array_push($sections,$row);
		}
		else if($prev_section == $current_section)
		{
			//remove last element and add current
			array_pop($sections);
			array_push($sections,$row);
		}
		else if($prev_section > $current_section)
		{
			$sections = array_slice($sections,0,$length);
			array_push($sections, $row);
		}
		else 
		{
			trigger_error("ERROR", E_USER_ERROR);
		}
	}
	
	//gene entry
	if(is_numeric($row[0]))
	{
		if($prev_section > $current_section)
		{
			//remove number of elements and add current
			$diff = $prev_section-$current_section;
			$sections = array_slice($sections,0,count($sections)-$diff);
		}
		
		$has_elements = true;
		
		$key = implode("->", $sections);
		if(!isset($data[$key]))	$data[$key] = array();
		$data[$key][] = $row;
	}
}

//convert data to TSV
$out_file = new Matrix();
$out_file->setHeaders(array("gene","details","source"));
$genes = array();
foreach($data as $section => $gs)
{
	foreach($gs as $g)
	{
		list($tmp) = explode(";", $g);
		list(,$tmp) = explode(" ", $tmp);
		
		$index = array_search($tmp,$genes);
		if($index!==false)
		{
			//duplicate entry
			$tmp_row = $out_file->get($index,1).":".$section;
			$out_file->set($index,1,$tmp_row);
			continue;
		}

		$genes[] = $tmp;
		$row = array($tmp,$section,"kegg");
		$out_file->addRow($row);
	}
}
$out_file->toTSV($out);