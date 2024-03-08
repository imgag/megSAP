<?php


error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("an_filter_dragen_somatic", "Add additional filter annotation to Dragen vcf.");
$parser->addInfile("in",  "VCF.GZ file generated from somatic Dragen calling", false);
$parser->addString("tumor_name",  "Tumor name", false);
$parser->addString("normal_name",  "Normal name", false);
$parser->addOutfile("out",  "VCF.GZ file with filters added in the filter column", false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

$genome = genome_fasta($build);


$vcf = $parser->tempFile("_unpacked.vcf");
$parser->exec("bgzip", "-d -c $in > $vcf", true);

//left-align
$vcf_aligned = $parser->tempFile("_aligned.vcf");
$parser->exec(get_path("ngs-bits")."VcfLeftNormalize", "-stream -in $vcf -out $vcf_aligned -ref $genome", true);

//sort
$vcf_sorted = $parser->tempFile("_sorted.vcf");
$parser->exec(get_path("ngs-bits")."VcfSort","-in $vcf_aligned -out $vcf_sorted", true);


//################################################################################################
//Filter variants
//################################################################################################

$variants = Matrix::fromTSV($vcf_sorted);
$variants_filtered = new Matrix();

//get indicies
$colnames = $variants->getHeaders();
$colidx_tumor = array_search($tumor_name, $colnames);
$colidx_normal = array_search($normal_name, $colnames);

//quality cutoffs (taken from Strelka)
$min_td = 20;
$min_taf = 0.05;
$min_tsupp = 3;
$min_nd = 20;
$max_naf_rel = 1/6;
$min_call_sq = 5;
$min_filter_sq = 17.5;

//set comments and column names
$filter_format = '#FILTER=<ID=%s,Description="%s">';
$comments = [
	sprintf($filter_format, "all-unknown", "Allele unknown"),
	sprintf($filter_format, "special-chromosome", "Special chromosome"),
	sprintf($filter_format, "depth-tum", "Sequencing depth in tumor is too low (< {$min_td})"),
	sprintf($filter_format, "freq-tum", "Allele frequency in tumor < {$min_taf}"),
	sprintf($filter_format, "depth-nor", "Sequencing depth in normal is too low (< {$min_nd})"),
	sprintf($filter_format, "freq-nor", "Allele frequency in normal > ".number_format($max_naf_rel, 2)." * allele frequency in tumor"),
	sprintf($filter_format, "lt-3-reads", "Less than {$min_tsupp} supporting tumor reads"),
	sprintf($filter_format, "weak-evidence", "Somatic Quality lower than {$min_filter_sq}")
	];

$variants_filtered->setComments(array_merge($variants->getComments(), $comments));
$variants_filtered->setHeaders($colnames);

$count_passing = 0;
for($i = 0; $i < $variants->rows(); ++$i)
{
	$row = $variants->getRow($i);

	$ref = $row[3];
	$alt = $row[4];
	$format = $row[8];
	$tumor = $row[$colidx_tumor];
	$normal = $row[$colidx_normal];

	$filters = [];

	$filter = array_diff(explode(";", $row[6]), ["."]);

	if (!preg_match("/^[acgtACGT]*$/", $alt))
	{
		$filter[] = "all-unknown";
	}
	if (chr_check($row[0], 22, false) === FALSE)
	{
		$filter[] = "special-chromosome";
	}
	$calls = [];
	
	//Somatic quality
	$f_parts = explode(":", $format);
	$i_sq_value = array_search("SQ", $f_parts);
	
	$parts = explode(":", $tumor);
	if($i_sq_value !== false)
	{
		$snp_q = floatval($parts[$i_sq_value]);
	}
	else
	{
		$snp_q = 0;
	}
	
	if ($snp_q < $min_call_sq) continue;
	
	if ($snp_q < $min_filter_sq) $filter[] = "weak-evidence";

	list($td, $tf) = vcf_dragen_var($format, $tumor, $alt);
	list($nd, $nf) = vcf_dragen_var($format, $normal, $alt);
	$calls[] = [ $alt, $td, $tf, $nd, $nf, $filter ];

	foreach ($calls as $call)
	{
		$variant = $row;
		list($alt, $td, $tf, $nd, $nf, $filter) = $call;
		$variant[4] = $alt;

		if ($td * $tf < $min_tsupp) $filter[] = "lt-3-reads";
		if ($td < $min_td) $filter[] = "depth-tum";
		if ($nd < $min_nd) $filter[] = "depth-nor";
		if ($tf < $min_taf) $filter[] = "freq-tum";
		if ($nf > $max_naf_rel * $tf && $nd*$nf >= 2) $filter[] = "freq-nor";

		if (empty($filter))
		{
			$filter[] = "PASS";
			$count_passing++;
		}
	}
	
	$variant[6] = implode(";", $filter);
	$variants_filtered->addRow($variant);
}
$vcf_filtered = $parser->tempFile("_filtered.vcf");
$variants_filtered->toTSV($vcf_filtered);

//remove invalid variant
$vcf_invalid = $parser->tempFile("_filtered_invalid.vcf");
$parser->exec(get_path("ngs-bits")."VcfFilter", "-remove_invalid -ref $genome -in $vcf_filtered -out $vcf_invalid", true);

//left-align
$vcf_aligned = $parser->tempFile("_aligned.vcf");
$parser->exec(get_path("ngs-bits")."VcfLeftNormalize", "-stream -in $vcf_invalid -out $vcf_aligned -ref $genome", true);
$final = $vcf_aligned;


//flag off-target variants
if (!empty($target))
{
	$vcf_offtarget = $parser->tempFile("_filtered.vcf");
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $final -mark off-target -reg $target -out $vcf_offtarget", true);
	$final = $vcf_offtarget;
}

//remove artefacts specific for processing system (blacklist)
//artefacts are caused e.g. by hairpin sequences when using enzymatic digestion
//how artefacts are determined is documented in /mnt/storage3/users/ahsturm1/Sandbox/2023_01_31_twist_indel_artefacts/
if (!empty($target))
{
	$artefact_vcf = repository_basedir()."/data/misc/enzymatic_digestion_artefacts/".basename($target, ".bed").".vcf";
	if (file_exists($artefact_vcf))
	{
		$vcf_with_artefacts = dirname($out)."/".basename($out, ".vcf.gz")."_with_enzymatic_artefacts.vcf.gz";
		$parser->exec("bgzip", "-c $final > $vcf_with_artefacts", true);
		$vcf_no_artefacts = $parser->tempFile("_filtered_no_artefacts.vcf");
		$parser->exec(get_path("ngs-bits")."VcfSubstract", "-in $final -in2 $artefact_vcf -out $vcf_no_artefacts");
		$final = $vcf_no_artefacts;
	}
}

//zip and index output file
$parser->exec("bgzip", "-c $final > $out", true);
$parser->exec("tabix", "-p vcf $out", true);


?>