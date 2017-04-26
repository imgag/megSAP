<?php
	/** 
		@page annotate_max_ent_loss
	*/

	require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
	error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);
	
	
	function get_reference_seq($ss_pos, $ss_chrom, $left_portion, $right_portion)
	{
		return get_ref_seq($ss_chrom,$ss_pos-$left_portion+1,$ss_pos+$right_portion);
	}
	
	function write_MaxEnt_input_file($oufile,$seq)
	{
		file_put_contents($oufile,"\n".$seq);
	}
	
	function calculate_MaxEnt_score($ss_seq,$calculate5ss,$parser)
	{
		$tmp_max_ent_in = $parser->tempFile("_maxEntin.txt");
		$tmp_max_ent_out = $parser->tempFile("_maxEntout.txt");
		file_put_contents($tmp_max_ent_in,"\n".$ss_seq);
		$script= ($calculate5ss) ? "score5.pl" : "score3.pl";
		$tmp_maxEnt_output_file="/mnt/users/ahlenzf1/php/test/Tmp_out.txt"; //TODO absolut path!
		shell_exec("cd /mnt/share/opt/alamut-batch-standalone-1.5.1/MaxEnt/ && perl $script ".realpath("$tmp_max_ent_in").">".$tmp_max_ent_out);
		$result_file= file($tmp_max_ent_out);
		$result_parts=explode("\t",$result_file[0]);
		return trim($result_parts[1]);
	}
	
	function combine_ref_mut_scores($maxEnt_ref_scores,$maxEnt_mut_scores)
	{
		return trim($maxEnt_ref_score)."\t".$maxEnt_mut_score;;
	}
	
	function exchange_base($reference_seq,$new_base,$pos)
	{
		if (($pos>=strlen($reference_seq))||($pos<0)) return $reference_seq;
		return substr_replace($reference_seq, $new_base, $pos, 1);
	}
	
	function get_obs_base($line)
	{
		list(,,,,$obs_base)=explode("\t",$line);
		return $obs_base;
	}
	
	function calculated_max_ent_seq($chrom,$pos,$dist_mut_to_splice_site,$obs_base,$antisense)
	{
		$exon_portion5ss = 3; //to detect a 5'-spice-site, 3 bases of exon are needed
		$pos_of_mut_5ss= 2-$dist_mut_to_splice_site;
		$intron_portion5ss = 6; //to detect a 5'-spice-site, 6 bases of the intron are needed (it is assumed the mutated base itself counts towards intron)
		$exon_portion3ss = 3; //to detect a 3'-spice-site, 3 bases of the exon are needed (to have an example, it is assumed the mutated base itself counts towards intron)
		$pos_of_mut_3ss= 19-$dist_mut_to_splice_site;
		$intron_portion3ss = 20; //to detect a 3'-spice-site, 19 bases of the intron are needed (it is assumed the mutated base counts itself towards intron)
		if ($antisense)
		{
			$ref_5ss_seq=get_reference_seq($pos+$dist_mut_to_splice_site,$chrom, $intron_portion5ss, $exon_portion5ss)."\n";	
			$ref_3ss_seq=get_reference_seq($pos+$dist_mut_to_splice_site,$chrom, $exon_portion3ss, $intron_portion3ss)."\n";		
			$mut_5ss_seq=exchange_base($ref_5ss_seq,$obs_base,$pos_of_mut_5ss);
			$mut_3ss_seq=exchange_base($ref_3ss_seq,$obs_base,$pos_of_mut_3ss);
			
			return array(
			"ref_5ss_seq" => rev_comp($ref_5ss_seq),
			"ref_3ss_seq" => rev_comp($ref_3ss_seq),
			"mut_5ss_seq" => rev_comp($mut_5ss_seq),
			"mut_3ss_seq" => rev_comp($mut_3ss_seq)
			);
		}
		else
		{
			$ref_5ss_seq=get_reference_seq($pos+$dist_mut_to_splice_site,$chrom, $exon_portion5ss, $intron_portion5ss)."\n";	
			$ref_3ss_seq=get_reference_seq($pos+$dist_mut_to_splice_site,$chrom, $intron_portion3ss, $exon_portion3ss)."\n";		
			$mut_5ss_seq=exchange_base($ref_5ss_seq,$obs_base,$pos_of_mut_5ss);
			$mut_3ss_seq=exchange_base($ref_3ss_seq,$obs_base,$pos_of_mut_3ss);
			
			return array(
			"ref_5ss_seq" => $ref_5ss_seq,
			"ref_3ss_seq" => $ref_3ss_seq,
			"mut_5ss_seq" => $mut_5ss_seq,
			"mut_3ss_seq" => $mut_3ss_seq
			);
		}
	}
	
	function ss_info_header_lines()
	{
		$ss_info_header_lines=array();
		$ss_info_header_lines[]="##INFO=<ID=distNearestSS,Number=1,Type=Integer,Description=\"Distance to nearest splice site.\">\n";
		$ss_info_header_lines[]="##INFO=<ID=nearestSSType,Number=1,Type=String,Description=\"Nearest splice site type (5'/3').\">\n";
		$ss_info_header_lines[]="##INFO=<ID=wtMaxEntScore,Number=1,Type=Float,Description=\"wild type sequence MaxEntScan score at nearest splice site.\">\n";
		$ss_info_header_lines[]="##INFO=<ID=varMaxEntScore,Number=1,Type=Float,Description=\"variant sequence MaxEntScan score at nearest splice site.\">\n";
		return $ss_info_header_lines;
	}
	
	function add_info_fields($line,$info_fields)
	{
		$line_parts=explode("\t",$line);
		$line_parts[7].=";".implode(";",$info_fields);
		return implode("\t",$line_parts);
	}
	
	function iterate_infile($in,$parser)
	{
		$outlines=array();	
		
		foreach(file($in) as $line)
		{
			if ((trim($line)=="")||(starts_with($line, "##")))//if meta info or whitespace only
			{
				
				//save unchanged line
				$outlines[]=$line;
			}
			else
			{
				if (starts_with($line, "#CHROM\t"))//if header line
				{
					$outlines=array_merge($outlines,ss_info_header_lines());
					$outlines[]=$line;
					continue;
				}
				
				
				
				$wtMaxEntScore="wtMaxEntScore=";
				$varMaxEntScore="varMaxEntScore=";
				
				list($chr,$start,)=explode("\t",$line);
				list($dist_mut_to_splice_site,$antisense,$three_prime_end)=find_nearest_exon($chr,$start,$parser);
				$distNearestSS="distNearestSS=$dist_mut_to_splice_site";
				if ($three_prime_end) $nearestSSType="nearestSSType=3'";
				else $nearestSSType="nearestSSType=5'";
				
				$obs_base=get_obs_base($line);
				$splice_site_seqs=calculated_max_ent_seq($chr,$start,$dist_mut_to_splice_site,$obs_base,$antisense);
				$outline=$line;
				if ($three_prime_end)
				{
					$wtMaxEntScore="wtMaxEntScore=".calculate_MaxEnt_score($splice_site_seqs["ref_3ss_seq"],false,$parser);
					$varMaxEntScore="varMaxEntScore=".calculate_MaxEnt_score($splice_site_seqs["mut_3ss_seq"],false,$parser);
				}
				else
				{
					$wtMaxEntScore="wtMaxEntScore=".calculate_MaxEnt_score($splice_site_seqs["ref_5ss_seq"],true,$parser);
					$varMaxEntScore="varMaxEntScore=".calculate_MaxEnt_score($splice_site_seqs["mut_5ss_seq"],true,$parser);
				}
				$outlines[]=add_info_fields($line,array($distNearestSS,$nearestSSType,$wtMaxEntScore,$varMaxEntScore));				
			}			
		}
		return $outlines;// array($header_lines,$ref_3ss_seqs,$ref_5ss_seqs,$mut_5ss_seqs,$mut_3ss_seqs);//hash
	}
	
	function find_nearest_exon($chr,$pos,$parser)
	{
		$cidx = get_path("ngs-bits")."/Cidx";
		$ucsc_exons = get_path("data_folder")."/dbs/UCSC/exons.bed";
		$range_start=$pos-100;//it is extremely unlikely for a mutation to influence a splice site farer away
		$range_end=$pos+100;//it is extremely unlikely for a mutation to influence a splice site farer away
		
		/*
		chr1:95699787-95699870  chr1    95699786        95699871        RWDD3_Exon_Nr1
		chr1:95699787-95699870  chr1    95699786        95699871        TMEM56-RWDD3_Exon_Nr1*/

		list($exons) = $parser->exec($cidx, "-in $ucsc_exons -pos $chr:$range_start-$range_end", false);
		$parsed_exons=array();

		foreach ($exons as $exon)
		{
			list(,,$start,$end,,$strand_symbol)=explode("\t", $exon);
			$parsed_exons[]=array($start,$end,$strand_symbol);
		}
		
		$antisense=true;
		$three_prime_end=true;
		if (count($parsed_exons)==0)
		{
			return 1000;//no existing splice site nearby
		}
		else 
		{
			//find the nearest one compared to input
			$min_dist=1000;
			foreach ($parsed_exons as $parsed_exon)
			{
				if ((abs($pos-$parsed_exon[0]))<abs($min_dist))
				{
					$nearest_splice_site=$parsed_exon;
					$min_dist=$parsed_exon[0]-$pos;
					$antisense=($parsed_exon[2]=="-");
					$three_prime_end=!$antisense;					
				}
				if ((abs($pos-$parsed_exon[1]))<abs($min_dist))
				{
					$nearest_splice_site=$parsed_exon;
					$min_dist=$parsed_exon[1]-$pos;
					$antisense=($parsed_exon[2]=="-");
					$three_prime_end=$antisense;
				}
			}
		}
		return array($min_dist,$antisense,$three_prime_end);
	}
	
	// parse command line arguments
	$parser = new ToolBase("annotate_max_ent_loss", "Annotates VCF file with MaxEnt prediction scores.");
	$parser->addInfile("in",  "Input file in VCF format.", false);
	$parser->addOutFile("out", "Annotated VCF Outfile.", false);
	extract($parser->parse($argv));
	
	$outlines=iterate_infile($in,$parser);
	file_put_contents($out,$outlines);

?>