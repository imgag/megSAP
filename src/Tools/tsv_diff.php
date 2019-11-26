<?php
/** 
	@page tsv_diff
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("tsv_diff", "Writes out mismatching lines of two TSVs.");
$parser->addInfile("in1", "First input TSV file.", false);
$parser->addInfile("in2", "Sedond input TSV file.", false);
$parser->addInt("keep_cols", "Always show the first n columns.", true, 3);
$parser->addOutfile("out", "Output file.", true);
extract($parser->parse($argv));

//add suffix to 
if (is_dir($in2) && file_exists($in2."/".$in1))
{
	$in2 = $in2."/".$in1;
}

$c_match = 0;
$c_mismatch = 0;
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
		print $line1;
		continue;
	}
	
	//header > print
	if (starts_with($line1, "#") && starts_with($line2, "#"))
	{
		print $line1;
		continue;
	}
	
	//header/content > compare
	if ($line1!=$line2)
	{
		$parts1 = explode("\t", $line1);
		$parts2 = explode("\t", $line2);
		
		//always print first n columns
		for ($j=0; $j<$keep_cols; ++$j)
		{
			if ($parts1[$j]!=$parts2[$j])
			{
				trigger_error("Column with index {$j} (keep_cols) is different in files!\n", E_USER_ERROR); 
			}
			
			if ($j!=0) print "\t";
			print $parts1[$j];
		}
		
		for ($j=$keep_cols; $j<count($parts1); ++$j)
		{
			$res = "";
			$entry1 = trim($parts1[$j]);
			$entry2 = trim($parts2[$j]);
			if ($entry1!=$entry2)
			{
				$res = trim($parts1[$j])." >> ".trim($parts2[$j]); 
			}
			print "\t$res";
		}
		print "\n";
		++$c_mismatch;
	}
	else
	{
		++$c_match;
	}
}
print "##matching lines: $c_match\n";
print "##mismatching lines: $c_mismatch\n";

?>