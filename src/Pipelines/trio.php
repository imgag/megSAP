<?php

/**
	@page trio
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function determine_index($name, $parts)
{
	$index = array_search($name, $parts);
	
	if ($index===FALSE || $index==-1)
	{
		trigger_error("Could not determine index of column '$name' in header line: ".implode(" ", $parts), E_USER_ERROR);
	}
	
	return $index;
}

//parse command line arguments
$parser = new ToolBase("trio", "Trio analysis pipeline.");
$parser->addInfile("f", "BAM file of father.", false, true);
$parser->addInfile("m", "BAM file of mother.", false, true);
$parser->addInfile("c", "BAM file of child (index).", false, true);
$parser->addString("out_folder", "Output folder name.", false);
//optional
$parser->addInfile("system",  "Processing system INI file used for all samples (created from NGSD via processed sample 'c' by default).", true);
$steps_all = array("vc", "an", "cn");
$parser->addString("steps", "Comma-separated list of steps to perform:\nvc=variant calling, an=annotation, cn=copy-number analysis.", true, implode(",", $steps_all));
$parser->addFlag("no_check", "Skip gender check of parents and parent-child correlation check (otherwise done before variant calling)");
$parser->addInt("min_dp", "Minimum depth in all three samples.", true, 10);
$parser->addFloat("max_af", "Maximum allele frequency in 1000G, ExAC and gnomAD database.", true, 0.01);
$parser->addInt("max_ngsd", "Maximum occurances in NGSD with same genotype as child.", true, 20);

extract($parser->parse($argv));

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//extract processing system information from DB
$sys = load_system($system, basename($c, ".bam"));
$target_file = $sys['type']!="WGS" ? $sys['target_file'] : get_path("data_folder")."/enrichment/ssHAEv6_2017_01_05.bed";
if ($target_file=="")
{
	trigger_error("Cannot perform trio analysis without target region (processing systems of {$c} is '".$sys["name_short"]."')!", E_USER_ERROR);
}

//prepare multi-sample paramters
$args_multisample = [
	"-bams $f $m $c",
	"-status control control affected",
	"-out_folder $out_folder",
	"-system $system",
	"-prefix trio",
	];

//variant calling
if (in_array("vc", $steps))
{
	//check parent-child correlation
	if (!$no_check)
	{
		$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in $f $c -mode bam -max_snps 4000", true);
		$correlation = explode("\t", $output[0][1])[3];
		if ($correlation<0.5)
		{
			trigger_error("The genotype correlation of father and child is {$correlation}; it should be above 0.5!", E_USER_ERROR);
		}
		$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in $m $c -mode bam -max_snps 4000", true);
		$correlation = explode("\t", $output[0][1])[3];
		if ($correlation<0.5)
		{
			trigger_error("The genotype correlation of mother and child is {$correlation}; it should be above 0.5!", E_USER_ERROR);
		}
	}
	
	//check gender of parents
	if (!$no_check)
	{
		list($stdout) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $f", true);
		$gender = trim(substr($stdout[count($stdout)-1], 7));
		if (starts_with($gender,"unknown"))
		{
			$parser->log("Could not check gender of father. It could not be determine from BAM file.");
		}
		else if (!starts_with($gender,"male"))
		{
			trigger_error("Gender of father is not male: '$gender'!", E_USER_ERROR);
		}
		
		list($stdout) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $m", true);
		$gender = trim(substr($stdout[count($stdout)-1], 7));
		if (starts_with($gender,"unknown"))
		{
			$parser->log("Could not check gender of mother. It could not be determine from BAM file.");
		}
		else if (!starts_with($gender,"female"))
		{
			trigger_error("Gender of mother is not female: '$gender'!", E_USER_ERROR);
		}
	}

	//variant calling with multi-sample pipeline
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps vc", true);	
}

//annotation
if (in_array("an", $steps))
{
	//determine gender of child
	$gender = "n/a";
	list($stdout) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $c", true);
	foreach($stdout as $line)
	{
		if (starts_with($line, "gender:"))
		{
			list(, $gender) = explode(":", $line);
		}
	}
	$child_is_male = trim($gender)=="male";
	
	//annotation with multi-sample pipeline
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps an", true);	
	
	//parse GSvar file - exract interesting variants and statistics
	$vars_high_depth = 0;
	$vars_mendelian_error = 0;
	$vars_rare = 0;
	$vars_denovo = 0;
	$vars_recessive = 0;
	$vars_hemizygous = 0;
	$vars_comphet = 0;
	$vars_hemizygous_chrx = 0;
		
	$sample_c = basename($c, ".bam");
	$sample_f = basename($f, ".bam");
	$sample_m = basename($m, ".bam");
	
	$annotations = array();
	$genes_comp_mother = array();
	$genes_comp_father = array();
	$gsvar = "$out_folder/trio.GSvar";
	$h = fopen($gsvar, "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="") continue;
		
		//skip comments
		if (starts_with($line, "##")) continue;
		
		//header columns (determine column indices)		
		if (starts_with($line, "#"))
		{
			if ($sample_c=="" || $sample_f=="" || $sample_m=="")
			{
				trigger_error("At least one of the sample names could not be determined! C={$sample_c} F={$sample_f} M={$sample_m}", E_USER_ERROR);
			}
			
			$parts = explode("\t", $line);
			
			$i_c = determine_index($sample_c, $parts);
			$i_f = determine_index($sample_f, $parts);
			$i_m = determine_index($sample_m, $parts);
			$i_quality = determine_index("quality", $parts);
			$i_1000g = determine_index("1000g", $parts);
			$i_exac = determine_index("ExAC", $parts);
			$i_gnomad = determine_index("gnomAD", $parts);
			$i_ngsd_hom = determine_index("ihdb_allsys_hom", $parts);
			$i_ngsd_het = determine_index("ihdb_allsys_het", $parts);
			$i_class = determine_index("classification", $parts);
			$i_filter = determine_index("filter", $parts);
			$i_genes = determine_index("gene", $parts);
			
			continue;
		}
		
		//parse content lines
		$parts = explode("\t", $line);
		$chr = $parts[0];
		$chr_is_autosome = $chr!="chrX" && $chr!="chrY" && $chr!="chrMT";
		$geno_c = $parts[$i_c];
		$geno_f = $parts[$i_f];
		$geno_m = $parts[$i_m];
		
		//count mendelian errors
		if (!contains($parts[$i_quality], "low_DP") && $chr_is_autosome)
		{
			++$vars_high_depth;
			
			//hom, hom => het/wt
			if ($geno_f=="hom" && $geno_m=="hom" && $geno_c!="hom") $vars_mendelian_error += 1;
			//hom, het => wt
			else if ($geno_f=="hom" && $geno_m=="het" && $geno_c=="wt") $vars_mendelian_error += 1;
			else if ($geno_f=="het" && $geno_m=="hom" && $geno_c=="wt") $vars_mendelian_error += 1;
			//het, wt  => hom
			else if ($geno_f=="het" && $geno_m=="wt" && $geno_c=="hom") $vars_mendelian_error += 1;
			else if ($geno_f=="wt" && $geno_m=="het" && $geno_c=="hom") $vars_mendelian_error += 1;
			//wt, wt  => het/hom
			else if ($geno_f=="wt" && $geno_m=="wt" && $geno_c!="wt") $vars_mendelian_error += 1;
		}
		
		//skip WT variants
		if ($geno_c=="wt") continue;
		
		//skip non-rare variants in public databases
		$skip = false;
		if ($parts[$i_1000g]>$max_af || $parts[$i_exac]>$max_af || $parts[$i_gnomad]>$max_af) $skip = true;
		
		//skip variant with too high NGSD counts
		if ($geno_c=="hom" && $parts[$i_ngsd_hom]>$max_ngsd) $skip = true;
		if ($geno_c=="het" && $parts[$i_ngsd_het]>$max_ngsd) $skip = true;
		
		//skip low depth variant
		$quality = explode(";", $parts[$i_quality]);
		foreach($quality as $entry)
		{
			if (starts_with($entry, "DP="))
			{
				$depths = explode(",", substr($entry, 3));
				if (min($depths)<$min_dp) $skip = true;
			}
		}
		
		//make sure we consider class 3/4/5 variant
		$class = $parts[$i_class];
		if ($class=="1" || $class=="2") $skip = true;
		$class_ge_3 = $class=="3" || $class=="4" || $class=="5";
		
		//skip non-interesting variants
		if ($skip && !$class_ge_3) continue;
		++$vars_rare;
		
		$tag = $parts[0].":".$parts[1]." ".$parts[3].">".$parts[4];
		
		//determine  allele frquency of each sample
		foreach($quality as $entry)
		{
			if (starts_with($entry, "AF="))
			{
				$afs = explode(",", substr($entry, 3));
				$i_min = min($i_c,$i_f,$i_m);
				$af_c = $afs[$i_c-$i_min];
				$af_f = $afs[$i_f-$i_min];
				$af_m = $afs[$i_m-$i_min];
			}
		}
		
		//CASE 1: DENOVO
		if ($geno_c!="wt" && $geno_f=="wt" && $geno_m=="wt" && $af_f<0.05 && $af_m<0.05)
		{
			//print "DENOVO: $tag C=$geno_c($af_c) F=$geno_f($af_f) M=$geno_m($af_m) $class ".implode(",", $quality)."\n";
			$annotations[$tag][] = "trio_denovo";
			++$vars_denovo;
		}
		
		//CASE 2: RECESSIVE
		if (!$child_is_male || $chr_is_autosome) //skip chrX/Y for males
		{
			if ($geno_c=="hom" && $geno_f=="het" && $geno_m=="het")
			{
				//print "RECESSIVE: $tag C=$geno_c($af_c) F=$geno_f($af_f) M=$geno_m($af_m) $class ".implode(",", $quality)."\n";
				$annotations[$tag][] = "trio_recessive";
				++$vars_recessive;
			}
		}
		
		//CASE 3.1: HEMIZYGOUS
		if ($chr_is_autosome)
		{
			if (($geno_c=="hom" && $geno_f=="het" && $geno_m=="wt")
				||
				($geno_c=="hom" && $geno_f=="wt" && $geno_m=="het"))
			{
				//print "HEMIZYGOUS: $tag C=$geno_c($af_c) F=$geno_f($af_f) M=$geno_m($af_m) $class ".implode(",", $quality)."\n";
				$annotations[$tag][] = "trio_hemizygous";
				++$vars_hemizygous;
			}
		}
		
		//CASE 3.2: HEMIZYGOUS CHRX
		if ($chr=="chrX" && $child_is_male)
		{
			if ($geno_c=="hom" && $geno_f=="wt" && $geno_m=="het")
			{
				//print "HEMIZYGOUS_CHRX: $tag C=$geno_c($af_c) F=$geno_f($af_f) M=$geno_m($af_m) $class ".implode(",", $quality)."\n";
				$annotations[$tag][] = "trio_hemizygous_chrX";
				++$vars_hemizygous_chrx;
			}
		}
		
		//CASE 4: COMPOUND-HETEROZYGOUS (mark genes)
		if ($chr_is_autosome)
		{
			if (contains($line, ":HIGH:") || contains($line, ":MODERATE:") || contains($line, ":LOW:")) //filter for LOW/MODERATE/HIGH impact - otherwise we get too many UTR and intronic variants
			{
				
				if ($geno_c=="het" && $geno_f=="het" && $geno_m=="wt")
				{
					$genes = explode(",", $parts[$i_genes]);
					foreach($genes as $gene)
					{
						$genes_comp_father[$gene][] = $tag;
					}
				}
				
				if ($geno_c=="het" && $geno_f=="wt" && $geno_m=="het")
				{
					$genes = explode(",", $parts[$i_genes]);
					foreach($genes as $gene)
					{
						$genes_comp_mother[$gene][] = $tag;
					}
				}
			}
		}
	}
	
	//CASE 4: COMPOUND-HETEROZYGOUS (intersect genes and mark variants)
	$genes_inter = array_intersect(array_keys($genes_comp_father), array_keys($genes_comp_mother));
	foreach($genes_inter as $gene)
	{
		foreach($genes_comp_father[$gene] as $tag)
		{
			$annotations[$tag][] = "trio_comp_f";
			++$vars_comphet;
		}
		foreach($genes_comp_mother[$gene] as $tag)
		{
			$annotations[$tag][] = "trio_comp_m";
			++$vars_comphet;
		}
	}
	
	//write output file
	$tmp = $parser->tempFile(".GSvar");
	$h2 = fopen($tmp, "w");
	rewind($h);
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="") continue;
		
		//update analysis type
		if (starts_with($line, "##ANALYSISTYPE="))
		{
			$line = "##ANALYSISTYPE=GERMLINE_TRIO";
		}
		
		//add filter info lines
		if (starts_with($line, "#chr"))
		{
			fwrite($h2, "##FILTER=trio_denovo=Trio analyis: Variant is de-novo in child.\n");
			fwrite($h2, "##FILTER=trio_recessive=Trio analyis: Variant is recessively inherited from parents.\n");
			fwrite($h2, "##FILTER=trio_hemizygous=Trio analyis: Variant is hemizygous.\n");
			fwrite($h2, "##FILTER=trio_hemizygous_chrX=Trio analyis: Variant is hemizygous (for males on chrX).\n");
			fwrite($h2, "##FILTER=trio_comp_m=Trio analyis: Variant is compound-heteroygous inherited from mother.\n");
			fwrite($h2, "##FILTER=trio_comp_f=Trio analyis: Variant is compound-heteroygous inherited from father.\n");			
		}
		
		//add annotations
		if ($line[0]!="#")
		{
			$parts = explode("\t", $line);
			$tag = $parts[0].":".$parts[1]." ".$parts[3].">".$parts[4];
			if (isset($annotations[$tag]))
			{
				$filters = explode(";", $parts[$i_filter]);
				foreach($annotations[$tag] as $anno)
				{
					$filters[] = $anno;
				}
				$parts[$i_filter] = implode(";", $filters);
				$line = implode("\t", $parts);
			}
		}
		
		fwrite($h2, "$line\n");
	}
	fclose($h);
	fclose($h2);
	$parser->moveFile($tmp, $gsvar);
	
	print "Medelian errors: ".number_format(100.0*$vars_mendelian_error/$vars_high_depth, 2)."% (of {$vars_high_depth} high-depth autosomal variants)\n";
	print "Rare or class 3/4/5 variants: {$vars_rare} (min_dp, max_af, max_ngsd, classification>2)\n";
	print "\n";
	print "Denovo variants: {$vars_denovo}\n";
	print "Recessive variants: {$vars_recessive}\n";
	print "Compound-heterozygous variants: {$vars_comphet}\n";
	print "Hemizygous variants: {$vars_hemizygous}\n";
	print "Hemizygous variants (chrX of males) {$vars_hemizygous_chrx}\n";
}

//copy-number
if (in_array("an", $steps))
{
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps cn", true);	
}

?>
