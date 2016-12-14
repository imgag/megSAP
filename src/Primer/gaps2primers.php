<?php

/** 
	@page list2fasta
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("gaps2primers", "Converts text input format to FASTA.");
$parser->addInfile("in",  "Input gaps BED file.", false);
$parser->addInfile("primer",  "Input primer TXT file.", false);
$parser->addOutfile("out",  "Output TXT file.", false);
extract($parser->parse($argv));

/*
################################### PAIR ###################################
INPUT: BRCA1_2F	BRCA1_2R	GGACGTTGTCATTAGTTCTTTGG	CTTCCCTAGTATGTAAGGTC	chr17	332
BLAST: BRCA1_2F - chr17:41276255-41276277 perc_aligned:100.00 length:23/23 e-value:3e-04
BLAST: BRCA1_2R - chr17:41275946-41275965 perc_aligned:100.00 length:20/20 e-value:0.008
################################### PAIR ###################################
INPUT: BRCA1_3F	BRCA1_3R	GCTCAAAGTTGAACTTATTCAC	CAAAAGCTAATAATGGAGCCAC	chr17	190
BLAST: BRCA1_3F - chr17:41267870-41267891 perc_aligned:100.00 length:22/22 e-value:8e-04
BLAST: BRCA1_3R - chr17:41267702-41267723 perc_aligned:100.00 length:22/22 e-value:8e-04
################################### PAIR ###################################
INPUT: BRCA1_5F	BRCA1_5R	GCCTTTTGAGTATTCTTTCTAC	TCCTACTGTGGTTGCTTCC	chr17	193
BLAST: BRCA1_5F - chr17:41258597-41258618 perc_aligned:100.00 length:22/22 e-value:8e-04
BLAST: BRCA1_5R - chr17:41258426-41258444 perc_aligned:100.00 length:19/19 e-value:0.031
*/
//load primers
$file = file($primer);
$prim = array();
$blast_line_count = 0;
foreach($file as $line)
{
	$line = trim($line);
	if (contains($line, "## PAIR ##"))
	{
		$blast_line_count = 0;
	}
	if (starts_with($line, "BLAST:"))
	{
		$blast_line_count += 1;
		
		list(, $name, , $pos) = explode(" ", $line);
		list($chr, $start, $end) = explode("-", strtr($pos,":", "-"));
		
		if ($start>$end)
		{
			trigger_error("Start greater than end in chromosomal range '$chr:$start-$end'", E_USER_ERROR);
		}
		
		if ($blast_line_count==1)
		{
			$tmp = array($chr, $start, $end, $name);
		}
		else if ($blast_line_count==2)
		{
			//make sure the first primer is the first (in reference genome order)
			if ($tmp[1]>$start)
			{
				$prim[] = array(array($chr, $start, $end, $name), $tmp);
			}
			else
			{
				$prim[] = array($tmp, array($chr, $start, $end, $name));
			}
		}
		else
		{
			trigger_error("More than two BLAST lines for one PAIR in primer file {$primer}!", E_USER_ERROR);
		}
	}
}

/*
  PTEN	61	chr10	89725179-89725239
  BRCA1	48	chr17	41258513-41258560
  BRCA1	192	chr17	41244481-41244529, 41258513-41258560, 41246006-41246084, 41246235-41246250
*/
//process gaps file
$file = file($in);
$output = array();
foreach($file as $line)
{
	if (trim($line)=="" || $line[0]=="#") continue;
	
	list($chr, $start, $end, $gene) = explode("\t", $line);
	$start += 1; //converte from 0-based to 1-based coordinates
	if ($start>$end)
	{
		trigger_error("Start greater than end in gap '$chr:$start-$end'", E_USER_ERROR);
	}
	$gene = trim($gene); //remove newline
	if($gene=="") $gene = "n/a";
		
	$output[] = "### GAP $gene $chr:$start-$end ###";
	
	foreach($prim as $primers)
	{
		list($p1, $p2) = $primers;
		if ($p1[0]==$chr)
		{
			$primer_start = $p1[2] + 1;
			$primer_end = $p2[1] - 1;
			
			if (range_overlap($start, $end, $primer_start, $primer_end))
			{
				$output[] = "PRIMER: ".$p1[3]." + ".$p2[3]." (".$p1[0].":$primer_start-$primer_end ~ ".number_format(100*(min($end, $primer_end)-max($start, $primer_start)+1)/($end-$start+1),2)."% of gap)";
			}
		}
	}
}

file_put_contents($out, implode("\n", $output));

?>