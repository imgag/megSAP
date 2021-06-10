<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("hgvs2vcf", "Convert a list of genomic HGVS.g variants to VCF format.");
$parser->addInfile("in", "Input TXT file (one variant per line).", false);
$parser->addOutfile("out",  "Output VCF file.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//write header
$out_h = fopen2($out, "w");
fwrite($out_h, "##fileformat=VCFv4.0\n");
fwrite($out_h, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
fwrite($out_h, "##INFO=<ID=ORIG,Number=A,Type=String,Description=\"Original HGVS variant string\">\n");
fwrite($out_h, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample\n");

//parse intput and convert
$file = file($in);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	$parts = explode("\t", strtr($line, array("g."=>"", ">"=>"\t", "_"=>"\t", ":"=>"\t", "dup"=>"\t", "del"=>"\t", "ins"=>"\t")));
	//print_r($parts);
	if (contains($line, ">"))
	{
		list($chr, $pos, $obs) = $parts;
		$ref = substr($pos, -1);
		$pos = substr($pos, 0, -1);
	}
	else if (contains($line, "ins") && contains($line, "del"))
	{
		list($chr, $pos, $pos2) = $parts;
		$obs = end($parts);
		$ref = get_ref_seq($build, $chr, $pos, $pos2);
	}
	else if (contains($line, "ins"))
	{
		list($chr, $pos, , $obs) = $parts;
		$anker = get_ref_seq($build, $chr, $pos, $pos);
		$ref = $anker;
		$obs = $anker.$obs;
	}
	else if (contains($line, "del"))
	{
		if (count($parts)==3)
		{
			list($chr, $pos, $ref) = $parts;
		}
		else
		{
			list($chr, $pos, , $ref) = $parts;
		}
		$pos -= 1;
		$anker = get_ref_seq($build, $chr, $pos, $pos);
		$obs = $anker;
		$ref = $anker.$ref;
	}
	else if (contains($line, "dup"))
	{
		if (count($parts)==3)
		{
			list($chr, $pos, $ref) = $parts;
		}
		else
		{
			list($chr, $pos, , $ref) = $parts;
		}
		$pos -= 1;
		$anker = get_ref_seq($build, $chr, $pos, $pos);
		$obs = $anker.$ref.$ref;
		$ref = $anker.$ref;
	}
	else
	{
		trigger_error("Invalid line $line", E_USER_ERROR);
	}
	fwrite($out_h, "$chr\t$pos\t.\t$ref\t$obs\t30\tPASS\tORIG=$line\tGT\t0/1\n");
}
fclose($out_h);

?>