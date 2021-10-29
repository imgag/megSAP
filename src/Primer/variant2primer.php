<?php

/** 
	@page variant2primer
*/

//parse command line arguments
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("variant2primer", "Finds primers for given region list");
$parser->addInfile("in",  "Input region list.", false);
$parser->addOutfile("out",  "Output primer list.", false);
$parser->addInt("min_len", "Minimum product length.", true, 200);
$parser->addInt("max_len", "Maximum product length.", true, 300);
$parser->addInt("num_of_res", "Maximal Number of Primer Pairs.", true, 5);
$parser->addInt("min_primer_size", "Minimal primer length.", true, 18);
$parser->addInt("max_primer_size", "Maximal primer length.", true, 24);
$parser->addInt("opt_primer_size", "Optimal Number of Primer Pairs.", true, 18);
extract($parser->parse($argv));

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function convert_region_list_to_bed($region_list_file,$out_file,$size)
{
	file_put_contents($out_file,"");//clear result file
	//build proper bed file from region list by converting "chr: start-end"-format lines
	//resulting bed entries are of size $size with given region in the center
	$region_lines=file($region_list_file);
	//print($region_lines);
	foreach($region_lines as $region_line)
	{
		
		list($chr,$start,$end)= preg_split("/[\s,\t,:,-]+/",trim($region_line));//split by space,tab,"-" and ":"
		$reg_length=$end-$start;
		$end=$end+ceil(($size-$reg_length)/2);
		$start=$start-floor(($size-$reg_length)/2);
		file_put_contents($out_file,implode("\t",array($chr,$start,$end))."\n", FILE_APPEND);
	}
	
}

function build_primer3_input_files($fasta_file,$directory,$min_len,$max_len,$num_of_res,$region_lines,$min_primer,$max_primer,$opt_primer)
{
	//build input files for primer 3 based on a fasta file, one for each sequence
	$fasta_lines=file($fasta_file);
	$outfiles=array();
	for ($index = 0; $index <count($fasta_lines); $index++) 
	{	
		if (($index+1) % 2 === 0)//if sequence line
		{
			$seq_ID=($index+1)/2;
			$seq=trim($fasta_lines[$index]);
			list(,$start,$end)= preg_split("/[\s,\t,:,-]+/",trim($region_lines[$index/2]));
			
			$query_length=$start-$end;
			$query_start=(int)(strlen($seq)/2)-(int)(($end-$start)/2);//substract 2 because of possible rounding errors
			$out_file=$directory."/".$seq_ID;
			$outfiles[]=$out_file;
			file_put_contents($out_file,"");//clear out file
			file_put_contents($out_file,"SEQUENCE_ID=".$seq_ID."\n", FILE_APPEND);
			file_put_contents($out_file,"SEQUENCE_TEMPLATE=".$seq."\n", FILE_APPEND);
			file_put_contents($out_file,"SEQUENCE_TARGET=".$query_start.",".($end-$start)."\n", FILE_APPEND);
			file_put_contents($out_file,"SEQUENCE_EXCLUDED_REGION=".($query_start-35).",".(($end-$start)+70)."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_TASK=pick_detection_primers"."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_PICK_LEFT_PRIMER=1"."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_PICK_RIGHT_PRIMER=1"."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_OPT_SIZE=".$opt_primer."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_MAX_SIZE=".$max_primer."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_MIN_SIZE=".$min_primer."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_MAX_NS_ACCEPTED=0"."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_PRODUCT_SIZE_RANGE=".max(($end-$start),$min_len)."-".$max_len."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_NUM_RETURN=".$num_of_res."\n", FILE_APPEND);
			file_put_contents($out_file,"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=".get_path("primer3")."primer3_config/"."\n", FILE_APPEND);
			file_put_contents($out_file,"="."\n", FILE_APPEND);
		}
	}
	return $outfiles;
}

function print_primer3_info()
{	
	//print primer3 information
	exec(get_path("primer3")."primer3_core --version 2>&1", $output, $error);
	$version="";
	foreach($output as $line)
	{
		if ((strpos($line, "libprimer3 release") !== false) && preg_match("/\d+\.\d+\.\d+/", $line, $hits))
		{
			$version = $line;
			break;
		}
	}
	print "Primer3 info: ".$version."\n";
}

function print_bedToFasta_info()
{	
	//print BedToFasta information
	exec(get_path("ngs-bits")."BedToFasta --version 2>&1", $output, $error);
	$version="";
	foreach($output as $line)
	{
		
		if (strpos($line, "BedToFasta") !== false)
		{
			$version = trim($line);
			break;
		}
	}
}


function parse_primer3_result($primer3_result_file)
{
	//parse primer3 result to data structure vector (common_line,list(primerpair1_lines,primerpair2_lines,...)
	$result_lines=file($primer3_result_file);
	$result_vector=array();
	$result_vector[0]=array();
	$index=0;
	$primer_pairs_found=false;
	while(($index<count($result_lines)))
	{
		if(substr($result_lines[$index], 0, 13 ) == "PRIMER_PAIR_0")
		{
			$primer_pairs_found=true;
			break;
		}
		$result_vector[0][]=$result_lines[$index]; 
		$index++;
	}
	if ($primer_pairs_found)
	{
		$primer_pair=0;
		$result_vector[$primer_pair+1]=array();
		foreach(array_slice($result_lines, $index,-1) as $result_line)//every line but the common ones at the beginning and the last one
		{
			$next_primer_pair=((substr($result_line, 0, 12 ) == "PRIMER_PAIR_")&&(substr($result_line, 12, strlen($primer_pair))!=$primer_pair));
			if ($next_primer_pair)
			{
				$primer_pair++;
				
				$result_vector[$primer_pair+1]=array();
			}
			$result_vector[$primer_pair+1][]=$result_line;
		}
	}
	return($result_vector);
}

//convert input to bed file
$temp_bed=$parser->tempFile("","bed");
$temp_fasta=$parser->tempFile("","fasta");
convert_region_list_to_bed($in,$temp_bed,$max_len*2);

//convert bed file (location) to fasta (sequences)
$parser->exec(get_path("ngs-bits")."BedToFasta", "-in ".$temp_bed." -out ".$temp_fasta." -ref ".get_path("GRCh37_data_folder")."/genomes/GRCh37.fa", true);

//run primer3 and convert output
$region_lines=file($in);
$tempDirectory=$parser->tempFolder();
$primer3_input_files=build_primer3_input_files($temp_fasta,$tempDirectory,$min_len,$max_len,$num_of_res,$region_lines,$min_primer_size,$max_primer_size,$opt_primer_size);
$index=0;
//data structure to store results, a map of query region -> array, where each array has one subsubarray to store common information 
//and one additional subsubarray per returned primer pair
$primer3result_map=array();
foreach($primer3_input_files as $primer3_file)
{
	$parser->exec(get_path("primer3")."primer3_core", "<$primer3_file >> $tempDirectory/$index.result" , true);
	$primer3result_map[trim($region_lines[$index])]=parse_primer3_result("$tempDirectory/$index.result");
	$index++;
}

//serialize the final map so it can be written to a file which is passed to check.php
$serialized_results=serialize($primer3result_map);
file_put_contents($out,$serialized_results);

?>