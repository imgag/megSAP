<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/functions.php");

$fcount = 0;
$pcount = 0;

$group = $argv[1];

for ($i=2; $i<$argc; ++$i)
{
	$a = $argv[$i];
	
	//filter by tool group
	if ($group!="ALL")
	{
		$text = file_get_contents($a);
		if (!contains($text, "/$group/"))
		{
			continue;
		}
	}
	
	print "TEST ".strtr($a, array("./"=>"", ".php"=>""));
	
	$time_start = microtime(true);
	$output = array();
	exec("php $a", $output);
	
	print " (".number_format(microtime(true) - $time_start, 2)."s)\n";

	foreach ($output as $line)
	{
		if (contains($line, " PASSED"))
		{
			++$pcount;
		}
		if (contains($line, " FAILED"))
		{
			++$fcount;
			print "  $line\n";
		}
	}	
}

//output
print "\n";
print "Summary:\n";
print "  - passed: $pcount\n";
print "  - failed: $fcount\n";

?>
