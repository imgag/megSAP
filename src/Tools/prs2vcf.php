<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("prs2vcf", "Convert a PRS file from https://www.pgscatalog.org/ into a VCF file.");
$parser->addInfile("in", "Input TXT file with PGS info.", false);
$parser->addOutfile("out",  "Output VCF file.", false);
//optional 
$parser->addString("build", "The genome build to use.", true, "GRCh37");
extract($parser->parse($argv));

//init
$genome_fasta = genome_fasta($build);

//write comment section
$temp_file = $parser->tempFile(".vcf");
$out_h = fopen2($temp_file, "w");
fwrite($out_h, "##fileformat=VCFv4.2\n");
fwrite($out_h, "##fileDate=".date("Ymd")."\n");

//parse input and convert
$header_parsed = false;
$file = file($in);
foreach($file as $line)
{
	if(starts_with($line, "#"))
	{
		// parse comment part
		if(starts_with($line, "# PGS ID"))
		{
			$pgs_id = trim(explode("=", $line)[1]);
			fwrite($out_h, "##pgs_id=$pgs_id\n");
		}
		else if(starts_with($line, "# Reported Trait"))
		{
			$trait = trim(explode("=", $line)[1]);
			fwrite($out_h, "##trait=$trait\n");
		}
		else if(starts_with($line, "# Original Genome Build"))
		{			
			// check if provided build and build in PRS file match
			$build_prs = trim(explode("=", $line)[1]);
			if($build_prs != $build) 
			{
				trigger_error("Provided genome build \'$build\' does not match genome build of PRS file (\'$build_prs\')!", E_USER_ERROR);
			}
			fwrite($out_h, "##build=$build\n");
		}
		else if(starts_with($line, "# Number of Variants"))
		{
			$n_var = trim(explode("=", $line)[1]);
			fwrite($out_h, "##n_var=$n_var\n");
		}
		else if(starts_with($line, "# PGP ID"))
		{
			$pgp_id = trim(explode("=", $line)[1]);
			fwrite($out_h, "##pgp_id=$pgp_id\n");
		}
		else if(starts_with($line, "# Citation"))
		{
			$citation = trim(explode("=", $line)[1]);
			fwrite($out_h, "##citation=$citation\n");
		}
	}
	else
	{
		if(!$header_parsed)
		{
			// parse header line
			$column_headers = explode("\t", trim($line));
			$chr_idx = array_search("chr_name", $column_headers);
			if($chr_idx === false) trigger_error("Mandatory column 'chr_name' is missing in input file!", E_USER_ERROR);
			$pos_idx = array_search("chr_position", $column_headers);
			if($pos_idx === false) trigger_error("Mandatory column 'chr_position' is missing in input file!", E_USER_ERROR);
			$ref_idx = array_search("reference_allele", $column_headers);
			if($ref_idx === false) trigger_error("Mandatory column 'reference_allele' is missing in input file!", E_USER_ERROR);
			$alt_idx = array_search("effect_allele", $column_headers);
			if($alt_idx  === false) trigger_error("Mandatory column 'effect_allele' is missing in input file!", E_USER_ERROR);
			$weight_idx = array_search("effect_weight", $column_headers);
			if($weight_idx  === false) trigger_error("Mandatory column 'effect_weight' is missing in input file!", E_USER_ERROR);

			//write header line
			fwrite($out_h, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n");
			
			$header_parsed = true;
			continue;
		}

		// parse data line
		$data_line = explode("\t", trim($line));
		
		$chr = $data_line[$chr_idx];
		$pos = $data_line[$pos_idx];
		$ref = $data_line[$ref_idx];
		$alt = $data_line[$alt_idx];
		$weight = $data_line[$weight_idx];

		// add "chr"-prefix
		if (!starts_with(strtolower(trim($chr)), "chr"))
		{
			$chr = "chr".strtoupper($chr);
		}

		// write VCF line
		fwrite($out_h, "$chr\t$pos\t.\t$ref\t$alt\t.\t.\tWEIGHT=$weight\n");

	}
	
}
fclose($out_h);

//check if all required header items are parsed:
if(!isset($pgs_id)) trigger_error("PGS ID missing in PRS file!");
if(!isset($trait)) trigger_error("Reported Trait missing in PRS file!");
if(!isset($build_prs)) trigger_error("Original Genome Build missing in PRS file!");
if(!isset($n_var)) trigger_error("Number of Variants missing in PRS file!");
if(!isset($pgp_id)) trigger_error("PGP ID missing in PRS file!");
if(!isset($citation)) trigger_error("Citation missing in PRS file!");


//left-align VCF file
$normalize_out = $parser->tempFile("_leftNormalized.vcf");
$parser->exec(get_path("ngs-bits")."VcfLeftNormalize", "-ref $genome_fasta -in $temp_file -out $normalize_out");
$temp_file = $normalize_out;



//calculate PRS for all WGS samples of the NGSD and calculate distribution (percentiles)
$distribution_file = $parser->tempFile("_distribution.tsv");
$parser->execTool("Tools/calculate_prs_distribution.php", "-in $temp_file -out $distribution_file");

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
$vcf_content = file($temp_file);
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
if ($header_idx < 0) trigger_error("VCF header not found in VCF '$temp_file'!", E_USER_ERROR);
array_splice($vcf_content, $header_idx, 0, $dist_header);

// write annotated VCF to disk
file_put_contents($out, $vcf_content);


//verify output file
$parser->exec(get_path("ngs-bits")."VcfCheck", "-in $out -lines 0 -ref $genome_fasta", true);


?>