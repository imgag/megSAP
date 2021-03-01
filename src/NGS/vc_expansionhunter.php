<?php

/**
  @page vc_expansionhunter
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_expansionhunter", "Call repeat expansions with Expansion Hunter. Creates an VCF file.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addOutfile("out", "Output VCF file.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addString("pid", "Processed sample name (e.g. 'GS120001_01'). If unset BAM file name will be used", true);
extract($parser->parse($argv));

// use BAM file name as fallback if no processed sample name is provided
if(!isset($pid)) $pid = basename($in);


// prepare command
$expansion_hunter_binary = get_path("expansion_hunter");
$args = [];
$args[] = "--reads $in";
$args[] = "--reference ".genome_fasta($build);
$args[] = "--output-prefix ".dirname($out)."/".basename($out, ".vcf");
		
//determine correct variant catalog for the current sample
$variant_catalog = dirname(get_path("expansion_hunter"), 2)."/variant_catalog/".strtolower($build)."/variant_catalog.json";
if (!file_exists($variant_catalog))
{
	trigger_error("No varaint catalog for given genome build '$build' found in '".dirname(get_path("expansion_hunter"), 2)."/variant_catalog/"."'!", E_USER_ERROR);
}
$args[] = "--variant-catalog $variant_catalog";
		
// add gender info
if (db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
	$info = get_processed_sample_info($db, $pid, false);
	if (!is_null($info))
	{
		if ($info['gender'] != 'n/a') $args[] = "--sex ".$info['gender'];
	}
}

// run Expansion Hunter
$parser->exec($expansion_hunter_binary, implode(" ", $args));


// add ExpansionHunter version string to VCF
list($stdout, $stderr, $ec) = $parser->exec($expansion_hunter_binary, "-v");

// parse STDERR
list($date, $tool_version) = explode(",", $stdout[0]);
$tool_version = strtr($tool_version, array("[Starting " => "", "]" => ""));

// read VCF file into memory
$vcf_file_content = file($out);

// write back to disk and insert additional comments
$fh = fopen($out, 'w');
$header_written = false;
foreach($vcf_file_content as $line)
{
	fwrite($fh, $line);
	// add info after the first line (##fileformat=...)
	if (!$header_written && starts_with($line, "##"))
	{
		// add file/call info to file
		fwrite($fh, "##filedate=$date\n");
		fwrite($fh, "##source={$tool_version}\n");
		fwrite($fh, "##reference=".genome_fasta($build)."\n");
		$header_written = true;
	}
	
}
?>