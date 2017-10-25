<?php

/**
	@page igv_session
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("igv_session", "Create an IGV session file.");

//mandatory arguments
$parser->addInfileArray("in", "List of files used as tracks.", false, true);
$parser->addOutfile("out", "Output IGV session XML file.", false);

//optional arguments
$parser->addFlag("relative", "Make paths relative to output directory.");
$parser->addString("genome", "IGV genome build.", true, "1kg_v37");
$parser->addString("locus", "Locus/region to display.", true, "all");

//extract arguments
extract($parser->parse($argv));

$out_lines = [];

$out_lines[] = <<<EOT
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{$genome}" hasGeneTrack="true" hasSequenceTrack="true" locus="{$locus}" version="8">
	<Resources>
EOT;

foreach ($in as $resource)
{
	if ($relative)
	{
		$resource = relative_path(dirname($out), $resource);
	}

	$out_lines[] = <<<EOT

		<Resource path="{$resource}" />
EOT;
}

$out_lines[] = <<<EOT

	</Resources>
	<HiddenAttributes>
		<Attribute name="NAME"/>
		<Attribute name="DATA FILE"/>
		<Attribute name="DATA TYPE"/>
	</HiddenAttributes>
</Session>

EOT;

file_put_contents($out, $out_lines);
