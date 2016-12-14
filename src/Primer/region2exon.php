<?php
/** 
	@page region2exon
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("region2exon", "Looks up exons for given regions.");
$parser->addInfile("in",  "Input 1-row bed file.", false);
$parser->addOutfile("out",  "Input bed file.", false);

extract($parser->parse($argv));

$query_regions = file($in);
$out_data=array();
foreach ($query_regions as $query_region)
{
	list($chr1,$start1,$end1)=preg_split("/[\s,\t,:,-]+/",$query_region);

	$cidx = get_path("ngs-bits")."/Cidx";
	$ucsc_exons = get_path("data_folder")."/dbs/UCSC/exons.bed";

	/*
	chr1:95699787-95699870  chr1    95699786        95699871        RWDD3_Exon_Nr1
	chr1:95699787-95699870  chr1    95699786        95699871        TMEM56-RWDD3_Exon_Nr1*/

	list($exons) = $parser->exec($cidx, "-in $ucsc_exons -pos $chr1:$start1-$end1", false);
	$parsed_exons=array();

	foreach ($exons as $exon)
	{
		list(,,$start,$end,$gene_exom,$ori)=explode("\t", $exon);
		$parsed_exons[]=array($start,$end,$gene_exom,$ori);
	}

	if (count($parsed_exons)==0)
	{
		$out_data[]="Not within an exon";
	}
	else 
	{
		if (count($parsed_exons)>1)//if multiple exones and/or genes exists at region
		{
			//check whether there are multiple exons or just the same used by multiple genes
			$synonymous=true;
			$start=$parsed_exons[0][0];
			$end=$parsed_exons[0][1];
			foreach ($parsed_exons as $parsed_exon)
			{
				if (($parsed_exon[0]!=$start)||($parsed_exon[1]!=$end))//at least on exon is not synonymous 
				{
					$synonymous=false;
					{
						$out_data[]="Non-synonymous exons in region";
					}
					break;
				}
			}
		}	
		//write all exons to file
		foreach ($parsed_exons as $parsed_exon)
		{
			$out_data[]=$chr1."\t".$parsed_exon[0]."\t".$parsed_exon[1]."\t".$parsed_exon[2]."\t".$parsed_exon[3];
		}
	}
}
file_put_contents($out,implode("\n",$out_data));
?>