<?php
/** 
	@page kasp_check 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function error_probabilty($match_prob, $matches, $mismatches)
{
	//calcualte value
	$p = 0.0;
	$count = $matches + $mismatches;
	for ($i=$matches; $i<=$count; ++$i)
	{
		$q = pow(1.0-$match_prob, $count-$i) * pow($match_prob, $i) * factorial($count) / factorial($i) / factorial($count - $i);
		$p += $q;
	}
	
	return $p;
}

///Converts a LightCycler text export file to TSV format
function txt2geno_lightcycler($snps, $in_file, $out_file)
{
	$skipped = array();
	
	//parse data
	$parsed = array();
	$file = file($in_file);
	foreach($file as $line)
	{
		//parse
		if (trim($line)=="" || starts_with($line, "Experiment") || starts_with($line, "Include")) continue;
		$parts = explode("\t", nl_trim($line));
		if (count($parts)<9)
		{
			trigger_error("Found line with less than 9 tab-separated parts '$in_file'!", E_USER_WARNING);
			continue;
		}
		list($include, $color, $position, $sample, $x, $y, $call, $score, $status) = $parts;
		
		//filter
		if (!isset($snps[$position[0]])) continue;
		if (!preg_match("/^[0-9]{6,}/", $sample) && !preg_match("/^DNA-[0-9]{6,}/", $sample) && !starts_with($sample, "FO"))
		{
			$skipped[$sample] = true;
			continue;
		}
		
		//determine genotype (KASP)
		list($rs, $chr, $pos, $ref, $alt, $reverse) = $snps[$position[0]];
		if (is_null($reverse))
		{
			$geno = "n/a";
		}
		else if ($call == "Allele X")
		{
			$geno = $reverse ? $alt : $ref;
		}
		else if ($call == "Allele Y")
		{
			$geno = $reverse ? $ref : $alt;
		}
		else if ($call == "Both Alleles")
		{
			$geno = $ref."/".$alt;
		}
		else
		{
			$geno = "n/a";
		}
		
		if (isset($parsed[$sample]) && isset($parsed[$sample][$rs]))
		{
			trigger_error("Found two entries for sample '$sample' and SNP '$rs' in '$in_file'!", E_USER_WARNING);
		}
		
		$parsed[$sample][$rs] = array("geno"=>$geno, "score"=>$score);
	}

	//write header
	$output = array();
	$outline = "#sample";
	foreach($snps as $row => $data)
	{
		list($rs, $chr, $pos, $x, $y) = $data;
		$outline .= "\t{$rs}_{$chr}_{$pos}_{$x}_{$y}";
	}
	$output[] = $outline;

	//write sample lines
	foreach($parsed as $sample => $sample_data)
	{
		$outline = $sample;
		foreach($snps as $row => $row_data)
		{
			$rs = $row_data[0];
			if (isset($sample_data[$rs]))
			{
				$outline .= "\t".$sample_data[$rs]["geno"];
			}
			else
			{
				$outline .= "\tn/a";
			}
		}
		$output[] = $outline;
	}
	
	//print skipped samples
	foreach($skipped as $sample => $dummy)
	{
		print "Skipped sample '{$sample}'. Only samples with DNA and FO number are processed!\n";
	}
	
	file_put_contents($out_file, implode("\n", $output));
}

///Converts a StepOnePlus text export file to TSV format
function txt2geno_steponeplus($snps, $in_file, $out_file)
{
	//parse data
	$parsed = array();
	$file = file($in_file);
	foreach($file as $line)
	{
		$line = trim($line);
		
		//only content lines
		if (!preg_match("/^[A-H][0-9]+\t/", $line)) continue;
		
		list($well, $sample, $rs, , , , , $quality, $call) = explode("\t", $line);
		list(, $rs) = explode("_", $rs);
		if (trim($call)=="" || starts_with($call, "Negative Control") || starts_with($call, "Undetermined"))
		{
			$call = "n/a";
		}
		else
		{
			list(, $call) = explode(" ", $call);
		}
		if (trim($call)=="") $call = "n/a";
		
		$parsed[$sample][$rs] = array("geno"=>$call, "score"=>$quality);
	}

	//write header
	$output = array();
	$outline = "#sample";
	foreach($snps as $row => $data)
	{
		list($rs, $chr, $pos, $x, $y) = $data;
		$outline .= "\t{$rs}_{$chr}_{$pos}_{$x}_{$y}";
	}
	$output[] = $outline;

	//write sample lines
	foreach($parsed as $sample => $sample_data)
	{
		$outline = $sample;
		foreach($snps as $row => $row_data)
		{
			list($rs, $chr, $pos, $ref, $alt, $reverse) = $row_data;
			if (!is_null($reverse) && isset($sample_data[$rs]))
			{
				$outline .= "\t".$sample_data[$rs]["geno"];
			}
			else
			{
				$outline .= "\tn/a";
			}
		}
		$output[] = $outline;
	}

	file_put_contents($out_file, implode("\n", $output));
}

//returns the genotype(s) for a sample at a certain position, or 'n/a' if the minimum depth was not reached.
function ngs_geno($bam, $chr, $pos, $ref, $min_depth)
{	
	//get pileup
	//print get_path("samtools")." mpileup -aa -f ".genome_fasta("GRCh38")." -r $chr:$pos-$pos $bam\n";
	list($output) = exec2(get_path("samtools")." mpileup -aa -f ".genome_fasta("GRCh38")." -r $chr:$pos-$pos $bam");
	//print_r($output);
	list($chr2, $pos2, $ref2, , $bases) = explode("\t", $output[0]);;
	
	//count bases
	$bases = strtoupper($bases);
	$counts = array("A"=>0, "C"=>0, "G"=>0, "T"=>0);
	for($i=0; $i<strlen($bases); ++$i)
	{
		$char = $bases[$i];
		if ($char=="." || $char==",") $char = $ref;
		if (isset($counts[$char]))
		{
			++$counts[$char];
		}
	}
	arsort($counts);
	
	//check depth
	$keys = array_keys($counts);
	$c1 = $counts[$keys[0]];
	$c2 = $counts[$keys[1]];
	if ($c1+$c2<$min_depth) return "n/a";
	
	//determine genotype
	$b1 = $keys[0];
	$b2 = ($c2>3 && $c2/($c1+$c2)>0.1) ? $keys[1] : $keys[0];
	
	return "$b1/$b2";
}


//searches for sample in NGSD and returns processed sample meta data
function sample_from_ngsd(&$db, $dna_number, $irp, $itp, $ibad)
{
	$output = array();
	
	$project_conditions = "(p.type='diagnostic'".($irp ? " OR p.type='research'" : "").($itp ? " OR p.type='test'" : "").")";
	
	//### get sample name ###

	if (!starts_with($dna_number, "FO"))
	{
		//1. try (DNA prefix (new schema))
		$res = $db->executeQuery("SELECT s.name FROM sample s WHERE s.name LIKE 'DNA{$dna_number}%' AND EXISTS (SELECT ps.id FROM processed_sample ps, project p WHERE ps.project_id=p.id AND ps.sample_id=s.id AND {$project_conditions})");
		
		//2. try (DX prefix and external name (old schema))
		if (count($res)==0)
		{
			$res = $db->executeQuery("SELECT s.name FROM sample s WHERE (s.name LIKE 'DX{$dna_number}%' OR s.name_external LIKE '%{$dna_number}%') AND EXISTS (SELECT ps.id FROM processed_sample ps, project p WHERE ps.project_id=p.id AND ps.sample_id=s.id AND {$project_conditions})");
		}
	}
	else //search for FO number
	{
		$res = $db->executeQuery("SELECT s.name FROM sample s WHERE (s.name LIKE '{$dna_number}%' OR s.name_external LIKE '%{$dna_number}%') AND EXISTS (SELECT ps.id FROM processed_sample ps, project p WHERE ps.project_id=p.id AND ps.sample_id=s.id AND {$project_conditions})");
		
		
		if (count($res)==0)
		{
			//2. try (without '-')
			$dna_number = str_replace("-", "", $dna_number);
			$res = $db->executeQuery("SELECT s.name FROM sample s WHERE (s.name LIKE 'DX{$dna_number}%' OR s.name_external LIKE '%{$dna_number}%') AND EXISTS (SELECT ps.id FROM processed_sample ps, project p WHERE ps.project_id=p.id AND ps.sample_id=s.id AND {$project_conditions})");
		}
	}
	
	//determine processed sample meta data
	$ngsbits = get_path("ngs-bits");
	foreach($res as $row)
	{
		$sample = $row['name'];
		list($stdout) = exec2("{$ngsbits}NGSDExportSamples -sample {$sample} ".($ibad ? "" : "-no_bad_samples")." -run_finished -add_path SAMPLE_FOLDER | {$ngsbits}TsvSlice -cols 'name,project_type,project_name,path,quality'");
		foreach($stdout as $line)
		{
			$line = trim($line);
			if ($line=="" || $line[0]=="#") continue;
			list($ps, $project_type, $project_name, $path, $quality) = explode("\t", $line);
			if ($project_type=="research" && !$irp) continue;
			if ($project_type=="test" && !$itp) continue;
			if ($project_name=="RPGR-Ex15") continue;
			$output[$sample][] = array($ps, $path, $quality);
		}
	}
	
	return $output;
}

//imports the result into NGSD
function import_ngsd(&$db, $ps, $rmp, $c_both, $c_match, $user_id)
{
	//determine processed sample id from BAM name
	$ps_id = get_processed_sample_id($db, $ps);
	
	$db->executeStmt("DELETE FROM kasp_status WHERE processed_sample_id=:ps_id", array("ps_id"=>$ps_id));
	$db->executeStmt("INSERT INTO kasp_status VALUES (:ps_id,:rmp,:c_both,:c_match,now(),:calculated_by)", array("ps_id"=>$ps_id, "rmp"=>$rmp, "c_both"=>$c_both, "c_match"=>$c_match, "calculated_by"=>$user_id));
}

//returns if a valid error probability exists in NGSD for the sample
function kasp_exists(&$db, $ps)
{
	$ps_id = get_processed_sample_id($db, $ps);
	$error_prob = $db->getValue("SELECT random_error_prob FROM kasp_status WHERE processed_sample_id={$ps_id}", -1);
	return $error_prob>=0 && $error_prob<=1;
}

//parse command line arguments
$parser = new ToolBase("kasp_check", "Checks that the genotype of a KASP essay and the NGS sample match.");
$parser->addInfile("in", "KASP text export file.", false);
$parser->addOutfile("out", "Output file in TSV format.", false);
$parser->addEnum("snps", "SNP set.", false, array("set1", "set2"));
$parser->addEnum("format", "Input text format.", false, array("LightCycler", "StepOnePlus"));
$parser->addInt("min_depth", "Minimal depth for genotype determination in NGS data", true, 10);
$parser->addInt("mep", "Maximum error probabilty of of genotype matches.", true, 0.01);
$parser->addFlag("pmm", "Prints mismatches to the command line.");
$parser->addFlag("irp", "Include research projects.");
$parser->addFlag("itp", "Include test projects.");
$parser->addFlag("ibad", "Include bad processed samples.");
$parser->addFlag("missing_only", "Only perform KASP<>NGS check for samples that don't have a entry in NGSD yet.");
$parser->addString("user", "Name of the user performing the KASP analysis and import (current user if unset).", true, "");
extract($parser->parse($argv));

//get user ID from NGSD
if ($user=="") $user = exec('whoami');
$db = DB::getInstance("NGSD");
$user_id = $db->getValue("SELECT id FROM user WHERE user_id='".$user."' AND active='1'", -1);
if ($user_id==-1)
{
	trigger_error("User '$user' not found in NGSD!", E_USER_ERROR);
}
print "$user $user_id\n";

if ($snps=="set1")
{
	$snps = array(
		"A" => array("rs6666954", "chr1", "78112493", "T", "C", true),
		"B" => array("rs4411641", "chr2", "146839405", "A", "G", false),
		"C" => array("rs11130795", "chr3", "60912762", "T", "C", false),
		"D" => array("rs6841061", "chr4", "185078389", "G", "A", false),
		"E" => array("rs37535", "chr5", "58321576", "G", "C", false),
		"F" => array("rs9388856", "chr6", "130827723", "A", "T", false),
		"G" => array("rs1393978", "chr8", "106224052", "G", "T", false),
		"H" => array("rs12682834", "chr9", "87447908", "A", "G", null), //null=skipped (this SNP produces quite a lot false calls)
		"I" => array("rs2583136", "chr11", "13081377", "G", "A", false),
		"J" => array("rs10748087", "chr12", "67801315", "C", "G", false),
		"K" => array("rs2988039", "chr13", "79192053", "A", "G", false),
		"L" => array("rs8045964", "chr16", "81783128", "C", "T", false),
		"M" => array("rs6074704", "chr20", "14186637", "A", "G", false),
		"N" => array("rs6512586", "chr20", "49684609", "G", "A", true),
	);
}
if ($snps=="set2")
{
	$snps = array(
		"A" => array("rs2275276", "chr1", "45508256", "G", "A", false),
		"B" => array("rs2274064", "chr1", "183573252", "T", "C", false),
		"C" => array("rs6788448", "chr3", "193491389", "T", "C", false),
		"D" => array("rs1371932", "chr4", "67914681", "A", "G", false),
		"E" => array("rs2046402", "chr4", "84841232", "T", "C", false),
		"F" => array("rs2071303", "chr6", "26091108", "T", "C", false),
		"G" => array("rs2229384", "chr6", "49457808", "C", "T", false),
		"H" => array("rs7742431", "chr6", "78969860", "A", "G", false),
		"I" => array("rs41288423", "chr6", "100718219", "G", "A", false),
		"J" => array("rs10756457", "chr9", "13150532", "T", "C", false),
		"K" => array("rs4290270", "chr12", "72022455", "A", "T", null), //null=skipped (this SNP produces very many false calls)
		"L" => array("rs6313", "chr13", "46895805", "G", "A", false),
		"M" => array("rs5351", "chr13", "77901178", "T", "C", false),
		"N" => array("rs140679", "chr15", "27527530", "C", "T", false),
	);
}

//(re-)create TSV file
$tsv_file = substr($in, 0, -4)."_converted.tsv";
if ($format=="LightCycler")
{
	txt2geno_lightcycler($snps, $in, $tsv_file);
}
else
{
	txt2geno_steponeplus($snps, $in, $tsv_file);
}

$db = DB::getInstance("NGSD");
$output = array();
$output[] = "#sample_name\trandom_match_prob\tgenos_kasp\tgenos_both\tgenos_match\tbam\n";
$file = file($tsv_file);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	$genotypes = explode("\t", $line);
	$name = trim($genotypes[0]);
	//extract DNA number
	preg_match("/[0-9]{6,}/", $genotypes[0], $matches);
	if (count($matches) == 1)
	{
		$dna_number = $matches[0];
		$genotypes = array_slice($genotypes, 1);
	}
	else
	{
		//check FO number
		preg_match("/FO[-]{0,1}[0-9]{5}/", $genotypes[0], $matches);
		if (count($matches) == 1)
		{
			$dna_number = $matches[0];
			$genotypes = array_slice($genotypes, 1);
		}
		else
		{
			print "error - Invalid DNA number '".$genotypes[0]."' given!\n";
			$output[] = "##Invalid DNA number '".$genotypes[0]."'\n";
			continue;
		}
	}
	 
	//determine BAM file of sample
	print "KASP: {$name} (sample search: $dna_number)\n";
	$res = sample_from_ngsd($db, $dna_number, $irp, $itp, $ibad);
	if (count($res)==0)
	{
		print "  error - Could not find any processed samples of the DNA number '$dna_number' in NGSD!\n";
		$output[] = "##no NGSD entry for sample '$dna_number'\n";
	}
	else
	{
		foreach($res as $sample => $ps_data)
		{
			print "  NGSD sample: {$sample}\n";
			
			$bams_found = 0;
			foreach($ps_data as list($ps, $folder, $quality))
			{
				$bam = realpath("$folder/{$ps}.bam");
				if (!file_exists($bam)) //fallback to CRAM
				{
					$bam = realpath("$folder/{$ps}.cram");
				}
				if (file_exists($bam))
				{
					++$bams_found;
					
					$messages = array();
					
					//check if missing
					if ($missing_only && kasp_exists($db, $ps))
					{
						print "  $ps - skipped because KASP is already in NGSD\n";
						continue;
					}
					
					//find genotype matches
					$c_kasp = 0;
					$c_both = 0;
					$c_match = 0;
					for($i=0; $i<count($genotypes); ++$i)
					{
						$g = $genotypes[$i];
						if ($g=="n/a") continue;
						++$c_kasp;
						$snp = $snps[chr(65+$i)];
						//print "  SNP: ".$snp[1].":".$snp[2]." ".$snp[3].">".$snp[4]."\n";
						//print "  KASP: $g\n";
						$g2 = ngs_geno($bam, $snp[1], $snp[2], $snp[3], $min_depth);
						if ($g2=="n/a") continue;
						++$c_both;
						//print "  NGS: $g2\n";
						
						//normalize genotypes
						$g = explode("/", $g);
						sort($g);
						$g = implode("/", array_unique($g));
						$g_hom = (strlen($g)==1);
						
						$g2 = explode("/", $g2);
						sort($g2);
						$g2 = implode("/", array_unique($g2));
						$g2_hom = (strlen($g2)==1);
						
						if ($g==$g2)
						{
							++$c_match;
						}
						else if ($pmm) //debug
						{
							$messages[] = "MISMATCH - snp:".implode("|", $snp)." kasp:$g ngs:$g2";
						}
						else 
						{
							if ($g_hom && $g2_hom) //WT vs HOM
							{
								$messages[] = "WARNING - snp:".implode("|", $snp)." kasp:$g ngs:$g2";
							}
						}
					}

					//determine overall match
					print "    ".basename2($bam)." ".($quality=="bad" ? "(bad)" : "")." kasp:$c_kasp both:$c_both match:$c_match";
					if ($c_both<=6)
					{
						$messages[] = "ERROR - too few common SNPs";
					}
					else
					{
						$rmp = number_format(error_probabilty(0.5, $c_match, $c_both-$c_match), 6);
						print " random_match_probabilty:$rmp";
						if($rmp>=$mep)
						{
							$messages[] = "ERROR - random_match_probabilty too high";
						}
						$output[] = "$dna_number\t$rmp\t$c_kasp\t$c_both\t$c_match\t$bam\n";
					
						//NGSD import of results
						import_ngsd($db, $ps, $rmp, $c_both, $c_match, $user_id);
					}
					print "\n";
					foreach($messages as $message)
					{
						print "      $message\n";
					}
				}
			}
			
			//error 
			if($bams_found==0)
			{
				print "  error - Could not find any BAM file for sample '$dna_number'!\n";
				$output[] = "##no BAM file for sample '$dna_number'\n";
			}
		}
	}
}

//write output
file_put_contents($out, $output);


?>