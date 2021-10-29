<?php

/** 
	@page blast
*/

//parse command line arguments
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("gap_closed", "Looks for overlaps between ab1 file and gaps.");
$parser->addInfileArray("ab1",  "Input AB1 files.", false);
$parser->addInfile("gap",  "Input gap BED file.", false);
$parser->addInfile("primers",  "Input primer BED file.", true);
$parser->addOutfile("out",  "Output closed gap.", false);
$parser->addFlag("hq",  "Only use high quality bases to close gaps.");
$parser->addInt("qual_thresh",  "use threshold value for quality check.", true);
$parser->addInt("qual_window",  "use window value for quality check.", true);
extract($parser->parse($argv));

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function parse_gap_list($gap_list_file)
{
	$output = array();
	$lines = file($gap_list_file);
	foreach($lines as $line)
	{
		if (trim($line)=="" || $line[0]=="#") continue;
		
		list($chr, $start, $end) = explode("\t", $line);
		$start += 1; //convert from 0-based to 1-based coordinates
		if ($start>$end)
		{
			$tmp = $end;
			$end = $start;
			$start = $tmp;
		}
		
		if(!isset($output[$chr])) $output[$chr] = array();
		$output[$chr][] = array($start, $end);
	}
	return $output;
}

function parse_blast_results($blast_chrom_to_regions, $blast_result_file, $primer_position_f, $primer_position_r)
{
	$blast_result=file($blast_result_file);
	foreach($blast_result as $line)
	{
		if($line[0]!='#') //if not header line
		{
			list(,$chrom,,,,,,,$first,$second)=array_map('trim',explode("\t",$line));
			$start=min($first,$second);
			$end=max($first,$second);
			$parts=explode("\t",$primer_position_f);
			if ($primer_position_f&&($parts[0]==$chrom))
			{
				$primer_start=min($parts[1],$parts[2])-1;
				$primer_end=max($parts[1],$parts[2]+1);
				//primer start<=product start<=primer end
				if(($primer_start<=$start)&&($primer_end>$start))
				{
					$start=$primer_end+1;
				}
				//primer start<=product end<=primer end
				if(($primer_start<$end)&&($primer_end>=$end))
				{
					$end=$primer_start-1;
				}
				$parts=explode("\t",$primer_position_r);
				$primer_start=min($parts[1],$parts[2])-1;
				$primer_end=max($parts[1],$parts[2])+1;
				//primer start<=product start<=primer end
				if(($primer_start<=$start)&&($primer_end>$start))
				{
					$start=$primer_end+1;
				}
				//primer start<=product end<=primer end
				if(($primer_start<$end)&&($primer_end>=$end))
				{
					$end=$primer_start;
				}
				
			}
			if(array_key_exists($chrom,$blast_chrom_to_regions))
			{
				$blast_chrom_to_regions[$chrom][]=array($start,$end);
			}
			else
			{
				$blast_chrom_to_regions[$chrom]=array(array($start,$end));
			}
		}
	}
	return $blast_chrom_to_regions;
}

function get_overlaps_regions($region1, $region2)
{
		//returns the start and end of the overlap of the input regions
		$start=max($region1[0],$region2[0]);
		$end=min($region1[1],$region2[1]);
		return array($start,$end); //overlapping region, possibly =< 1 basepairs
}

function get_overlaps_region_lists($gaplist, $blastlist)
{
	//calculates the bases of each gap region which overlaps with at least one region of the blast results
	//returns dictionary of gap_string->numbers of overlapping bases
	$gap_closer=array();
	$gap_bitmap=array();
	foreach($gaplist as $gap)
	{
		$gap_string=implode('-',array($gap[0],$gap[1]));
		$size_of_gap=$gap[1]-$gap[0]+1;
		$gap_bitmap[$gap_string]=array_fill(0, $size_of_gap, false);//at beginning, no base is closed
		
		foreach($blastlist as $blastmatch)
		{
			$overlap=get_overlaps_regions($gap,$blastmatch);
			if ($overlap[0]<=$overlap[1])//if the overlap is not negative
			{
				$start_of_overlap_gap=$overlap[0]-$gap[0];//first base of gap which overlaps
				$end_of_overlap_gap=$overlap[1]-$gap[0];//first base of gap which overlaps
				for ($i=$start_of_overlap_gap;$i<=$end_of_overlap_gap;$i++)//set base positions covered to true
				{
					$gap_bitmap[$gap_string][$i]=true;
				}
			}
		}
		$gap_closer[$gap_string]=0;
		for ($jjj=0;$jjj<count($gap_bitmap[$gap_string]);$jjj++)//count number of covered base positions
		{
			if ($gap_bitmap[$gap_string][$jjj])
			{
				$gap_closer[$gap_string]++;
			}
		}
	}
	return $gap_closer;
}

function get_overlaps_result_dics($gap_results,$blast_results)
{
	//calculates the start and end of the longest overlap of each gap region of each chromosome with a region of the blast results
	//returns dictionary of chromomesomes->(gap_string->length of longest overlap)
	$gap_closer=array();
	foreach(array_keys($gap_results) as $chrom)
	{
		if(array_key_exists($chrom,$blast_results))
		{
			$gap_closer[$chrom]=get_overlaps_region_lists($gap_results[$chrom], $blast_results[$chrom]);
		}
		else
		{
			$gap_closer[$chrom]=get_overlaps_region_lists($gap_results[$chrom], array());//empty blast results for that chromosome
		}
	}
	return $gap_closer;
}

function add_gap_length($comparison_result)
{
	$closed_gaps=array();
	foreach($comparison_result as $chrom => $gaps)
	{
		foreach ($gaps as $gap=>$closed_length)
		//if ($closed_length>0)
		{
			list($gap_start,$gap_end)=explode("-",$gap);
			$gap_length=$gap_end-$gap_start+1;
			$closed_gaps[$chrom][$gap]=array($gap_length,$closed_length);
		}
	}
	return $closed_gaps;
}

//print TraceTuner information
$fasta_from_ab1_file=$parser->tempfile("from_ab1_file.FASTA");
$output=array();
exec(get_path("tracetuner")." 2>&1", $output, $error);
foreach($output as $line)
{
	if ((strpos($line, "Version") !== false) && preg_match("/\d+\.\d+\.\d+/", $line, $hits))
	{
		$version = $line;
		break;
	}
}
print "TraceTuner info: ".$version."\n";

//print Blast information
$output = $parser->exec("blastn", "-version", false);
print "BLAST info: ".$output[0][0]."\n";

$blast_chrom_to_regions=array();
$index=0;
if (!(isset($qual_window))) $qual_window=2;
if (!(isset($qual_thresh))) $qual_thresh=23;

if (isset($primers))//if a primer position file was provided
{
	$primer_positions=file($primers);
	if (count($primer_positions)!=2*count($ab1))
	{
		print "<b>ERROR! Number of Primers is not two times the number of ab1 files: ".count($ab1)." primers:".count($primer_positions)."</b>";
		exit;
	}
	
}
else
{
	$primer_positions=array();
}

$counter=0;
$primer_position_f="";
$primer_position_r="";
foreach ($ab1 as $ab1_file)
{
	//convert ab1 to FASTA
	$parser->exec(get_path("tracetuner"), "-nocall -trim_threshold $qual_thresh -trim_window $qual_window -sa ".$fasta_from_ab1_file." \"".$ab1_file."\"", true);

	//cut low quality region if flag is set
	if ($hq)
	{
		$fasta_file_lines=file($fasta_from_ab1_file);
		list($name,$num_of_bases,$start,$length)=explode(" ",trim($fasta_file_lines[0]));//get start and length of high quality region from FASTA header
		$fasta_seq=implode("",array_slice($fasta_file_lines,1));//get sequence from FASTA file
		$hq_fasta_seq=substr(preg_replace('/\s+/', '',$fasta_seq),$start,$length);//extract high quality sequence
		file_put_contents($fasta_from_ab1_file,implode(" ",array($name,$length,1,$length))."\n".$hq_fasta_seq);//write high quality fasta
	}
	//blast FASTA
	$blast_result_file=$parser->tempfile("blast_results.txt");
	$index++;
	$parser->exec("blastn", "-task blastn-short -db ".get_path("GRCh37_data_folder")."/dbs/blast/hg19 -query ".$fasta_from_ab1_file." -outfmt 7 -out ".$blast_result_file." -num_threads 8", true);

	//extract blast results, cut primer regions if available
	if ($primer_positions) 
	{
		$primer_position_f=$primer_positions[$counter];
		$primer_position_r=$primer_positions[$counter+1];
	}
	$blast_chrom_to_regions=parse_blast_results($blast_chrom_to_regions, $blast_result_file,$primer_position_f,$primer_position_r);
	
	$gaps_chrom_to_regions=parse_gap_list($gap);
	$overlaps=get_overlaps_result_dics($gaps_chrom_to_regions,$blast_chrom_to_regions);
	
	$closed_gaps=add_gap_length($overlaps);
	$counter+=2;
	
}
//output results
file_put_contents($out,"");//make sure output file is empty at beginning
foreach($closed_gaps as $chrom => $gaps)
{
	foreach($gaps as $gap => $bases_closed)	
	{
		$outstring=implode("\t",array($chrom, $gap, $bases_closed[0], $bases_closed[1]."\n"));
		file_put_contents($out,$outstring, FILE_APPEND);
	}
}
?>