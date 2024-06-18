<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("prs2vcf", "Convert a PRS file from https://www.pgscatalog.org/ into a VCF file.");
$parser->addOutfile("out",  "Output VCF file.", false);
//optional 
$parser->addInfile("pgs", "Input TXT file with PGS info.", true);
$parser->addInfile("vcf", "Input VCF file with PGS info.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFlag("skip_percentiles", "Skip percentile computation");
$parser->addString("exclude_disease_group", "Name of a disease group which should be excluded for percentile calculation", true, "");
$parser->addString("processing_system", "Processing system short name to which the percentile calculation should be limited to.", true, "");
$parser->addFlag("add_tsv", "Add an additional tsv output");
extract($parser->parse($argv));

//init
$genome_fasta = genome_fasta($build, ($build=="GRCh38"), ($build=="GRCh38")); //local data works only for GRCh38

if (!(isset($pgs) xor isset($vcf)))
{
	trigger_error("Tool needs either PGS file or VCF as input!", E_USER_ERROR);
}

if (isset($pgs))
{
	/*
		use PGS file as input
	*/

	//write comment section
	$temp_file = $parser->tempFile(".vcf");
	$out_h = fopen2($temp_file, "w");
	fwrite($out_h, "##fileformat=VCFv4.2\n");
	fwrite($out_h, "##fileDate=".date("Ymd")."\n");

	//parse input and convert
	$header_parsed = false;
	$file = file($pgs);
	foreach($file as $line)
	{
		if(starts_with($line, "#"))
		{
			// parse comment part
			if(starts_with($line, "# PGS ID") || starts_with($line, "#pgs_id="))
			{
				$pgs_id = trim(explode("=", $line)[1]);
				fwrite($out_h, "##pgs_id=$pgs_id\n");
			}
			else if(starts_with($line, "# Reported Trait") || starts_with($line, "#trait_reported="))
			{
				$trait = trim(explode("=", $line)[1]);
				fwrite($out_h, "##trait=$trait\n");
			}
			else if(starts_with($line, "# Original Genome Build") || starts_with($line, "#genome_build="))
			{			
				$build_prs = trim(explode("=", $line)[1]);
			}
			else if(starts_with($line, "# Number of Variants") || starts_with($line, "#variants_number="))
			{
				$n_var = trim(explode("=", $line)[1]);
				fwrite($out_h, "##n_var=$n_var\n");
			}
			else if(starts_with($line, "# PGP ID") || starts_with($line, "#pgp_id="))
			{
				$pgp_id = trim(explode("=", $line)[1]);
				fwrite($out_h, "##pgp_id=$pgp_id\n");
			}
			else if(starts_with($line, "# Citation") || starts_with($line, "#citation="))
			{
				$citation = trim(explode("=", $line)[1]);
				fwrite($out_h, "##citation=$citation\n");
			}
			else if(starts_with($line, "#HmPOS_build="))
			{			
				$build_prs_hm = trim(explode("=", $line)[1]);
			}
		}
		else
		{
			if(!$header_parsed)
			{
				// check if provided build and build in PRS file match
				$use_harmonized_coordinates = false;
				if($build_prs != $build) 
				{
					//PRS genome build doesn't match, check harmonized build
					if($build_prs_hm != $build) 
					{
						trigger_error("Provided genome build '$build' does not match genome build of PRS file (PRS build:'{$build_prs}', Harmonized PRS build: '{$build_prs_hm}')!", E_USER_ERROR);
					}
					else
					{
						//Harmonized build matches
						trigger_error("Provided genome build '$build' matches harmonized genome build of PRS file. Using harmonized coordinates.", E_USER_NOTICE);
						$use_harmonized_coordinates = true;
					}
				}
				fwrite($out_h, "##build=$build\n");

				// parse header line
				$column_headers = explode("\t", trim($line));
				if ($use_harmonized_coordinates)
				{
					$chr_idx = array_search("hm_chr", $column_headers);
					if($chr_idx === false) trigger_error("Mandatory column 'hm_chr' is missing in input file!", E_USER_ERROR);
					$pos_idx = array_search("hm_pos", $column_headers);
					if($pos_idx === false) trigger_error("Mandatory column 'hm_pos' is missing in input file!", E_USER_ERROR);
					//optional fields
					$rsid_idx = array_search("hm_rsID", $column_headers);
					if($rsid_idx === false) trigger_error("Optional column 'hm_rsID' is missing in input file. It will be missing in the output file.", E_USER_WARNING);
				}
				else
				{
					$chr_idx = array_search("chr_name", $column_headers);
					if($chr_idx === false) trigger_error("Mandatory column 'chr_name' is missing in input file!", E_USER_ERROR);
					$pos_idx = array_search("chr_position", $column_headers);
					if($pos_idx === false) trigger_error("Mandatory column 'chr_position' is missing in input file!", E_USER_ERROR);
					//optional fields
					$rsid_idx = array_search("rsID", $column_headers);
					if($rsid_idx === false) trigger_error("Optional column 'rsID' is missing in input file. It will be missing in the output file.", E_USER_WARNING);
				}
				
				$other_idx = array_search("other_allele", $column_headers);
				if($other_idx === false) 
				{
					$other_idx = array_search("reference_allele", $column_headers);
					if($other_idx === false) trigger_error("Mandatory column 'reference_allele' or 'other_allele' is missing in input file!", E_USER_ERROR);
				}
				$effect_idx = array_search("effect_allele", $column_headers);
				if($effect_idx  === false) trigger_error("Mandatory column 'effect_allele' is missing in input file!", E_USER_ERROR);
				$weight_idx = array_search("effect_weight", $column_headers);
				if($weight_idx  === false) trigger_error("Mandatory column 'effect_weight' is missing in input file!", E_USER_ERROR);
				$popaf_idx = array_search("allelefrequency_effect", $column_headers);
				if($popaf_idx  === false) trigger_error("Mandatory column 'allelefrequency_effect' is missing in input file!", E_USER_ERROR);

				//add info columns
				fwrite($out_h, "##INFO=<ID=WEIGHT,Number=1,Type=Float,Description=\"PRS weight of this variant.\">\n");
				fwrite($out_h, "##INFO=<ID=POP_AF,Number=1,Type=Float,Description=\"Population allele frequency of this variant.\">\n");
				fwrite($out_h, "##INFO=<ID=OTHER_ALLELE,Number=1,Type=String,Description=\"The other allele(s) at the loci. Note: this does not necessarily need to correspond to the reference allele.\">\n");
				fwrite($out_h, "##INFO=<ID=REF_IS_EFFECT_ALLELE,Number=0,Type=Flag,Description=\"Set if reference base is effect allele\">\n");
				//write header line
				fwrite($out_h, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n");
				
				$header_parsed = true;
				continue;
			}

			// parse data line
			$data_line = explode("\t", trim($line));
			
			$chr = $data_line[$chr_idx];
			// add "chr"-prefix
			if (!starts_with(strtolower(trim($chr)), "chr"))
			{
				$chr = "chr".strtoupper($chr);
			}
			$pos = $data_line[$pos_idx];
			$rsid = ($rsid_idx === false)?".":$data_line[$rsid_idx];
			$other = $data_line[$other_idx];
			$ref = get_ref_seq($build, $chr, $pos, intval($pos) + (strlen($other)-1), 0, ($build=="GRCh38"));
			$alt = $data_line[$effect_idx];

			//fix insertions:
			if((strlen($alt) > 1) && ($alt != $other))
			{
				//replace first base with ref
				$alt = $ref[0].substr($alt, 1);
			}
			$info_annotation = "";
			//wt variants
			if ($alt == $ref) 
			{
				$alt = $other; // use other as alt
				$info_annotation = ";REF_IS_EFFECT_ALLELE";
			}
			$weight = $data_line[$weight_idx];
			$popaf = $data_line[$popaf_idx];

			// write VCF line
			if ($other == $ref) fwrite($out_h, "$chr\t$pos\t.\t$ref\t$alt\t.\t.\tPOP_AF={$popaf};WEIGHT={$weight}{$info_annotation}\n");
			else fwrite($out_h, "$chr\t$pos\t.\t$ref\t$alt\t.\t.\tPOP_AF={$popaf};WEIGHT={$weight};OTHER_ALLELE={$other}{$info_annotation}\n");

		}
		
	}
	fclose($out_h);

	//check if all required header items are parsed:
	if(!isset($pgs_id)) trigger_error("PGS ID missing in PRS file!", E_USER_WARNING);
	if(!isset($trait)) trigger_error("Reported Trait missing in PRS file!", E_USER_WARNING);
	if(!isset($build_prs)) trigger_error("Original Genome Build missing in PRS file!", E_USER_WARNING);
	if(!isset($n_var)) trigger_error("Number of Variants missing in PRS file!", E_USER_WARNING);
	if(!isset($pgp_id)) trigger_error("PGP ID missing in PRS file!", E_USER_WARNING);
	if(!isset($citation)) trigger_error("Citation missing in PRS file!", E_USER_WARNING);

	// set input file for left normalization
	$input_vcf = $temp_file;
}
else
{
	/*
		use VCF file as input
	*/

	//parse header and check if all required meta data is available
	$file = file($vcf);
	foreach($file as $line)
	{
		if(starts_with($line, "##"))
		{
			if(starts_with($line, "##pgs_id="))
			{
				$pgs_id = trim(explode("=", $line)[1]);
				continue;
			}
			if(starts_with($line, "##trait="))
			{
				$trait = trim(explode("=", $line)[1]);
				continue;
			}
			if(starts_with($line, "##build="))
			{
				$build_prs = trim(explode("=", $line)[1]);
				if($build_prs != $build) 
				{
					trigger_error("Provided genome build \'$build\' does not match genome build of input VCF file (\'$build_prs\')!", E_USER_ERROR);
				}
				continue;
			}
			if(starts_with($line, "##n_var="))
			{
				$n_var = trim(explode("=", $line)[1]);
				continue;
			}
			if(starts_with($line, "##pgp_id="))
			{
				$pgp_id = trim(explode("=", $line)[1]);
				continue;
			}
			if(starts_with($line, "##citation="))
			{
				$citation = trim(explode("=", $line)[1]);
				continue;
			}
		}
		else
		{
			// VCF comment section parsed -> abort
		break;
		}
	}

	//check if all required header items are parsed:
	if(!isset($pgs_id)) trigger_error("PGS ID missing in PRS file!");
	if(!isset($trait)) trigger_error("Reported Trait missing in PRS file!");
	if(!isset($build_prs)) trigger_error("Original Genome Build missing in PRS file!");
	if(!isset($n_var)) trigger_error("Number of Variants missing in PRS file!");
	if(!isset($pgp_id)) trigger_error("PGP ID missing in PRS file!");
	if(!isset($citation)) trigger_error("Citation missing in PRS file!");

	// set input file for left normalization
	$input_vcf = $vcf;
}



//left-align VCF file
$annotate_gnomad_af = (isset($pgs) && ($build=="GRCh38"));
$normalize_out = $parser->tempFile("_leftNormalized.vcf");
$pipeline = array();
$pipeline[] = array(get_path("ngs-bits")."VcfLeftNormalize", "-stream -ref $genome_fasta -in $input_vcf");
$pipeline[] = array(get_path("ngs-bits")."VcfStreamSort", (($annotate_gnomad_af)?"":"-out {$normalize_out}"));
if ($annotate_gnomad_af)
{
	//gnomAD annotation
	$args_gnomad = array();
	$args_gnomad[] = "-out {$normalize_out}";
	$args_gnomad[] = "-source ".get_path("data_folder")."/dbs/gnomAD/gnomAD_genome_v3.1.2_GRCh38.vcf.gz";
	$args_gnomad[] = "-prefix gnomADg";
	$args_gnomad[] = "-info_keys AC,AF,Hom,Hemi,Het,Wt,AFR_AF,AMR_AF,EAS_AF,NFE_AF,SAS_AF";
	$args_gnomad[] = "-allow_missing_header";
	$pipeline[] = array(get_path("ngs-bits")."VcfAnnotateFromVcf", implode(" ", $args_gnomad));
} 
$parser->execPipeline($pipeline, "VCF normalize and sort");

if ($annotate_gnomad_af)
{
	//calculate diff between gnomad AF (EUR) and POP_AF
	$vcf_content = file($normalize_out);
	$vcf_out_content = array();
	$outliers = array();
	foreach ($vcf_content as $line) 
	{
		if (starts_with($line, "##"))
		{
			$vcf_out_content[] = $line;
			continue;
		}
		if (starts_with($line, "#"))
		{
			// add info comment before header
			$vcf_out_content[] = "##INFO=<ID=AF_DIFF,Number=1,Type=String,Description=\"Absolute difference between gnomadAF (gnomADg_NFE_AF) and population AF (POP_AF).\">\n";
			$vcf_out_content[] = $line;
			continue;
		}
		if (trim($line) == "")
		{
			//keep empty lines
			$vcf_out_content[] = $line;
			continue;
		}

		//calc diff on variant lines
		$columns = explode("\t", $line);
		$info_values = explode(";", trim($columns[7]));
		//extract af
		$gnomad_af = 0;
		$pop_af = 0;
		foreach ($info_values as $info_value) 
		{
			if (starts_with($info_value, "POP_AF="))
			{
				$pop_af = (float) explode("=", $info_value)[1];
				continue;
			}
			if (starts_with($info_value, "gnomADg_NFE_AF="))
			{
				$gnomad_af = (float) explode("=", $info_value)[1];
				continue;
			}
		}
		$diff = abs($gnomad_af-$pop_af);
		$info_values[] = "AF_DIFF=".$diff;

		if ($diff > 0.1) trigger_error("The given population AF of the PRS variant ".$columns[0].":".$columns[1]." ".$columns[3].">".$columns[4]." differs more than 0.1 ({$diff}) from the gnomAD AF!", E_USER_WARNING);

		$columns[7] = implode(";", $info_values);
		$vcf_out_content[] = implode("\t", $columns)."\n";
	}

	file_put_contents($normalize_out, $vcf_out_content);

}


$vcf_content = file($normalize_out);
if(!$skip_percentiles)
{
	//calculate PRS for all WGS samples of the NGSD and calculate distribution (percentiles)
	$distribution_file = $parser->tempFile("_distribution.tsv");
	$parser->execTool("Tools/calculate_prs_distribution.php", "-in $input_vcf -out $distribution_file".(($exclude_disease_group == "")?"":" -exclude_disease_group ".$exclude_disease_group).(($processing_system == "")?"":" -processing_system ".$processing_system));

	// parse distribution file
	$distribution = Matrix::fromTSV($distribution_file);

	// check if table is valid
	if ($distribution->cols() != 102)
	{
		trigger_error("PRS distribution file has to contain 102 columns, got ".$distribution->cols()."!", E_USER_ERROR);
	}
	if ($distribution->rows() < 1)
	{
		trigger_error("PRS distribution file does not contain any distribution!", E_USER_ERROR);
	}

	$pgs_idx = $distribution->getColumnIndex("pgs_id");
	$sample_count_idx = $distribution->getColumnIndex("sample_count");
	$distribution_found = false;
	$percentiles = array();
	$sample_count = 0;
	for ($row_idx=0; $row_idx < $distribution->rows(); $row_idx++) 
	{ 
		// skip entries which does not fit to the given PGS id
		if ($distribution->get($row_idx, $pgs_idx) != $pgs_id) continue;
		
		// parse line and exit
		$distribution_found = true;
		$sample_count = $distribution->get($row_idx, $sample_count_idx);
		$percentiles = array_slice($distribution->getRow($row_idx), 2);
		if (count($percentiles) != 100)
		{
			trigger_error("Only ".count($percentiles)." percentile values in line, expected 100!", E_USER_ERROR);
		}
		break;
	}

	$dist_header = array("##sample_count=$sample_count\n", "##percentiles=".implode(",", $percentiles)."\n");

	// add distribution to VCF
	//find header
	$header_idx = -1;
	for ($i=0; $i < count($vcf_content); $i++) 
	{ 
		if (starts_with($vcf_content[$i], "#CHROM	POS	ID	REF	ALT"))
		{
			$header_idx = $i;
		break;
		}
	}
	if ($header_idx < 0) trigger_error("VCF header not found in VCF '$input_vcf'!", E_USER_ERROR);
	array_splice($vcf_content, $header_idx, 0, $dist_header);
}
// write annotated VCF to disk
file_put_contents($out, $vcf_content);

if ($add_tsv)
{
	$tsv_out = substr($out,0,-3)."tsv";
	$parser->exec(get_path("ngs-bits")."VcfToTsv", "-in {$out} -out {$tsv_out}");
}

//verify output file
list($stdout, $stderr, $return) = $parser->exec(get_path("ngs-bits")."VcfCheck", "-in $out -lines 0 -ref $genome_fasta", true, true);
print "VcfCheck output: \n".implode("\n", $stdout)."\n";



?>