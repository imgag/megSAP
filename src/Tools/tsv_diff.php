<?php
/** 
	@page tsv_diff
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("tsv_diff", "Writes out mismatching lines of two TSVs.");
$parser->addInfile("in1", "First input TSV file.", false);
$parser->addInfile("in2", "Second input TSV file.", false);
$parser->addInt("keep_cols", "Always show the first n columns.", true, 3);
$parser->addOutfile("out", "Output file.", true);
extract($parser->parse($argv));

//add suffix to 
if (is_dir($in2) && file_exists($in2."/".$in1))
{
	$in2 = $in2."/".$in1;
}

if ($out=="") $out = 'php://stdout';
$h_o = fopen2($out, 'w');

$c_match = 0;
$c_mismatch = 0;
$c_changes_by_col = [];
$file1 = file($in1);
$file2 = file($in2);
for($i=0; $i<count($file1); ++$i)
{
	$line1 = $file1[$i];
	$line2 = $file2[$i];
	
	//empty line > skip
	if (trim($line1)=="" && trim($line2)=="")
	{
		continue;
	}
	
	//comment > print
	if (starts_with($line1, "##") && starts_with($line2, "##"))
	{
		fputs($h_o, $line1);
		continue;
	}
	
	//header > print
	if (starts_with($line1, "#") && starts_with($line2, "#"))
	{
		fputs($h_o, $line1);
		$headers = explode("\t", nl_trim(substr($line1, 1)));
		continue;
	}

	// check for header mismatches
	if (starts_with($line1, "#") || starts_with($line2, "#"))
	{
		trigger_error("Number of header lines is different in files (in line $i)!\n", E_USER_ERROR);
	}
	
	//header/content > compare
	$matching = true;
	if ($line1!=$line2)
	{
		$parts1 = explode("\t", $line1);
		$parts2 = explode("\t", $line2);
		
		//always print first n columns
		for ($j=0; $j<$keep_cols; ++$j)
		{
			if ($parts1[$j]!=$parts2[$j])
			{
				trigger_error("Column with index {$j} (keep_cols) is different in files (in line $i)!\n", E_USER_ERROR); 
			}
			
			if ($j!=0) fputs($h_o, "\t");
			fputs($h_o, $parts1[$j]);
		}
		
		for ($j=$keep_cols; $j<count($parts1); ++$j)
		{
			$res = "";
			$entry1 = trim($parts1[$j]);
			$entry2 = trim($parts2[$j]);
			
			//ignore order and whitespaces in some fields
			if ($entry1!=$entry2)
			{
				$col_name = $headers[$j];
				if ($col_name=="gene" || $col_name=="variant_type" || $col_name=="coding_and_splicing")
				{
					$entry1 = explode(",", $entry1);
					$entry1 = array_map('trim', $entry1);
					$entry1 = array_filter($entry1, function($str) { return $str!="";} );
					sort($entry1);
					$entry1 = implode(",", $entry1);
					
					$entry2 = explode(",", $entry2);
					$entry2 = array_map('trim', $entry2);
					$entry2 = array_filter($entry2, function($str) { return $str!="";} );
					sort($entry2);
					$entry2 = implode(",", $entry2);
				}
			}
			
			if ($entry1!=$entry2)
			{
				$res = $entry1." >> ".$entry2;
				$matching = false;
				@$c_changes_by_col[$j] += 1;
			}
			fputs($h_o, "\t$res");
		}
		fputs($h_o, "\n");
	}
	
	if ($matching)
	{
		++$c_match;
	}
	else
	{
		++$c_mismatch;
	}
}
print "matching lines: $c_match\n";
print "mismatching lines: $c_mismatch\n";
foreach($c_changes_by_col as $index => $count)
{
	print "changes in column '".$headers[$index]."': {$count}\n";
}

?>