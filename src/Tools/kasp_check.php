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
	//parse data
	$parsed = array();
	$file = file($in_file);
	foreach($file as $line)
	{
		//parse
		if (trim($line)=="" || starts_with($line, "Experiment") || starts_with($line, "Include")) continue;
		list($include, $color, $position, $sample, $x, $y, $call, $score, $status) = explode("\t", nl_trim($line));

		//filter
		if (!preg_match("/^[0-9]{6}/", $sample) && !preg_match("/^DNA-[0-9]{6}/", $sample) && !starts_with($sample, "FO")) continue;
		if (!isset($snps[$position[0]])) continue;
		
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

	file_put_contents($out_file, implode("\n", $output));
}

//returns the genotype(s) for a sample at a certain position, or 'n/a' if the minimum depth was not reached.
function ngs_geno($bam, $chr, $pos, $ref, $min_depth)
{	
	//get pileup
	list($output) = exec2(get_path("samtools")." mpileup -aa -r $chr:$pos-$pos $bam");
	list($chr2, $pos2, $ref2, , $bases) = explode("\t", $output[0]);;
	
	//count bases
	$bases = strtoupper($bases);
	$counts = array("A"=>0, "C"=>0, "G"=>0, "T"=>0);
	for($i=0; $i<strlen($bases); ++$i)
	{
		$char = $bases[$i];
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


//imports the result into NGSD
function sample_from_ngsd(&$db, $sample_name, $irp)
{
	//sample name
	$res = $db->executeQuery("SELECT s.name FROM sample s WHERE s.name='{$sample_name}' AND EXISTS (SELECT ps.id FROM processed_sample ps, project p WHERE ps.project_id=p.id AND ps.sample_id=s.id AND (p.type='diagnostic'".($irp ? " OR p.type='research'" : "")."))");
	if (count($res)>0) return $res;
	
	//DNA number
	if (!starts_with($sample_name, "FO"))
	{
		if (starts_with($sample_name, "DNA-"))
		{
			$sample_name = substr($sample_name, 4);
		}
		$sample_name = substr($sample_name, 0, 6);
		$res = $db->executeQuery("SELECT s.name FROM sample s WHERE (s.name='DX{$sample_name}' OR s.name_external LIKE '%{$sample_name}%') AND EXISTS (SELECT ps.id FROM processed_sample ps, project p WHERE ps.project_id=p.id AND ps.sample_id=s.id AND (p.type='diagnostic'".($irp ? " OR p.type='research'" : "")."))");
		if (count($res)>0) return $res;
	}
	
	//FO number
	$res = $db->executeQuery("SELECT s.name FROM sample s WHERE (s.name_external LIKE '%{$sample_name}%') AND EXISTS (SELECT ps.id FROM processed_sample ps, project p WHERE ps.project_id=p.id AND ps.sample_id=s.id AND (p.type='diagnostic'".($irp ? " OR p.type='research'" : "")."))");
	
	return $res;
}

//imports the result into NGSD
function import_ngsd(&$db, $bam, $rmp, $c_both, $c_match)
{
	//determine processed sample id from BAM name
	$ps_name = basename($bam,".bam");
	$ps_id = get_processed_sample_id($ps_name);
	
	$db->executeStmt("DELETE FROM kasp_status WHERE processed_sample_id=:ps_id", array("ps_id"=>$ps_id));
	$db->executeStmt("INSERT INTO kasp_status VALUES (:ps_id,:rmp,:c_both,:c_match)", array("ps_id"=>$ps_id, "rmp"=>$rmp, "c_both"=>$c_both, "c_match"=>$c_match));
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
extract($parser->parse($argv));

if ($snps=="set1")
{
	$snps = array(
		"A" => array("rs6666954", "chr1", "78578177", "T", "C", true),
		"B" => array("rs4411641", "chr2", "147596973", "A", "G", false),
		"C" => array("rs11130795", "chr3", "60898434", "T", "C", false),
		"D" => array("rs6841061", "chr4", "185999543", "G", "A", false),
		"E" => array("rs37535", "chr5", "57617403", "G", "C", false),
		"F" => array("rs9388856", "chr6", "131148863", "A", "T", false),
		"G" => array("rs1393978", "chr8", "107236280", "G", "T", false),
		"H" => array("rs12682834", "chr9", "90062823", "A", "G", null), //null=skipped (this SNP produces quite a lot false calls)
		"I" => array("rs2583136", "chr11", "13102924", "G", "A", false),
		"J" => array("rs10748087", "chr12", "68195095", "C", "G", false),
		"K" => array("rs2988039", "chr13", "79766188", "A", "G", false),
		"L" => array("rs8045964", "chr16", "81816733", "C", "T", false),
		"M" => array("rs6074704", "chr20", "14167283", "A", "G", false),
		"N" => array("rs6512586", "chr20", "48301146", "G", "A", true),
	);
}
if ($snps=="set2")
{
	$snps = array(
		"A" => array("rs2275276", "chr1", "45973928", "G", "A", false),
		"B" => array("rs2274064", "chr1", "183542387", "T", "C", false),
		"C" => array("rs6788448", "chr3", "193209178", "T", "C", false),
		"D" => array("rs1371932", "chr4", "68780399", "A", "G", false),
		"E" => array("rs2046402", "chr4", "85762385", "T", "C", false),
		"F" => array("rs2071303", "chr6", "26091336", "T", "C", false),
		"G" => array("rs2229384", "chr6", "49425521", "C", "T", false),
		"H" => array("rs7742431", "chr6", "79679577", "A", "G", false),
		"I" => array("rs41288423", "chr6", "101166095", "G", "A", false),
		"J" => array("rs10756457", "chr9", "13150531", "T", "C", false),
		"K" => array("rs4290270", "chr12", "72416235", "A", "T", null), //null=skipped (this SNP produces very many false calls)
		"L" => array("rs6313", "chr13", "47469940", "G", "A", false),
		"M" => array("rs5351", "chr13", "78475313", "T", "C", false),
		"N" => array("rs140679", "chr15", "27772676", "C", "T", false),
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
	$sample_name = trim($genotypes[0]);
	$genotypes = array_slice($genotypes, 1);
	
	//determine BAM file of sample
	$res = sample_from_ngsd($db, $sample_name, $irp);
	if (count($res)==0)
	{
		print "error - Could not find sample '$sample_name' in NGSD!\n";
	}
	foreach($res as $res_entry)
	{
		$sample = $res_entry["name"];
		print "$sample_name (".$res_entry["name"].")\n";
		
		//find BAM file(s) of sample
		$project_folder  = get_path("project_folder");
		list($stdout) = $parser->exec("find", "-L {$project_folder}/diagnostic/ ".($irp ? "{$project_folder}/research/" : "")." -maxdepth 3 -name \"{$sample}_0?.bam\" | grep -v bad | sort", true);	
		if (count($stdout)==0)
		{
			trigger_error("Could not find BAM file for sample $sample!", E_USER_ERROR);
		}	
		
		//check each matching BAM file
		foreach($stdout as $bam)
		{
			$messages = array();
			
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
			print "  $bam kasp:$c_kasp both:$c_both match:$c_match";
			if ($c_both<=6)
			{
				$messages[] = "ERROR - too few common snps";
				
				//handle case where KASP was ok, but sample is not
				if ($c_kasp>=10)
				{
					//NGSD import of results
					import_ngsd($db, $bam, 999, $c_both, $c_match);
				}
			}
			else
			{
				$rmp = number_format(error_probabilty(0.5, $c_match, $c_both-$c_match), 6);
				print " random_match_probabilty:$rmp";
				if($rmp>=$mep)
				{
					$messages[] = "ERROR - random_match_probabilty too high";
				}
				$output[] = "$sample_name\t$rmp\t$c_kasp\t$c_both\t$c_match\t$bam\n";
			
				//NGSD import of results
				import_ngsd($db, $bam, $rmp, $c_both, $c_match);
			}
			print "\n";
			foreach($messages as $message)
			{
				print "    $message\n";
			}
		}
	}
}

//write output
file_put_contents($out, $output);


?>