<?php
/** 
	@page primer2SNP
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("primer2SNP", "Looks up SNPs in a primer pair.");
$parser->addInfile("in",  "Input TXT file that contains two lines (space-separated chromosome, start position and end position).", false);
$parser->addOutfile("out",  "Output TXT file.", false);
//optional
$parser->addFloat("min_freq", "Minimum allele frequency of SNPs.", true, 0.0);
extract($parser->parse($argv));

//filters entry according to minimum frequency. Returns 'false' if frequency too low.
function filter_by_af($line, $min_freq)
{
	$parts = explode("\t", $line);
	if (count($parts)<8) return false;
	
	$af = 0.0;
	$info = explode(";", $parts[7]);
	foreach($info as $entry)
	{
		if (starts_with($entry, "AF="))
		{
			$af = substr($entry, 3); 
		}
	}
	if ($af<$min_freq) return false;
	$af = number_format(max(0.0001, $af), 4);
	
	$id = $parts[2];
	if ($id==".") $id = "";
	
	return "chr{$parts[0]}\t{$parts[1]}\t{$parts[1]}\t{$parts[3]}\t{$parts[4]}\t{$af}\t{$id}";
}

//init
$dbs = array(
	"dbSNP" => get_path("GRCh37_data_folder")."/dbs/1000G/1000g_v5b.vcf.gz",
	"ExAC" => get_path("GRCh37_data_folder")."/dbs/ExAC/ExAC_r0.3.1.vcf.gz",
	"gnomAD" => get_path("GRCh37_data_folder")."/dbs/gnomAD/gnomAD_genome_r2.1.1.vcf.gz"
);
$output = array();
$file = file($in);
$primers = array("forward" => $file[0], "reverse" => $file[1]);

//process
foreach($primers as $p_name => $p_data)
{
	list($chr, $start, $end) = explode(" ", $p_data);
	if (!starts_with($chr, "chr")) trigger_error("Chromosome '$chr' does not start with 'chr'", E_USER_ERROR);
	foreach($dbs as $db_name => $db_file)
	{
		list($snps) = $parser->exec("tabix", "$db_file ".substr($chr, 3).":$start-$end", false);
		foreach($snps as $snp)
		{
			$snp = filter_by_af($snp, $min_freq);
			if ($snp!==FALSE)
			{
				$output[] = "$p_name	$db_name	$snp";
			}
		}
	}
}

//output
file_put_contents($out, implode("\n", $output));

?>