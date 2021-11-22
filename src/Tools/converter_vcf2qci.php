<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_vcf2qci", "Generate QCI compliant files from strelka vcf files.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addString("t_id",  "Tumor ID in VCF file.", false);
$parser->addString("n_id",  "Normal ID in VCF file.", false);
$parser->addOutfile("out",  "Output file in VCF format.", false);
// optional
$parser->addFlag("pass", "Keep only variants that passed all filter criteria");
$parser->addString("build", "Build used for vcf file generation.", true, "GRCh38");
extract($parser->parse($argv));

$filec = Matrix::fromTSV($in);

$idx_tumor = $filec->getColumnIndex($t_id);
$idx_normal = $filec->getColumnIndex($n_id);
if(is_null($idx_tumor)) trigger_error("Could not identify tumor column in $in.", E_USER_ERROR);
if(is_null($idx_normal)) trigger_error("Could not identify normal column in $in.", E_USER_ERROR);

$file_qci = new Matrix();
$tmp_headers = $filec->getHeaders();
$tmp_headers[$idx_tumor] = "TUMOR";
$file_qci->setHeaders($tmp_headers);

// add additional format fields GT and AD
for($i=0; $i<$filec->rows();++$i)
{
	$row = $filec->getRow($i);	
	if($pass && $row[6]!="PASS" && $row[6] != "freq-tum")	continue; 
	
	//  format field: add allele frequency "AD - 0.3"
	$alleles = "0,0";
	
	$variant = "SNV";
	list(,,,$ref,$alt,,,,) = $row;
	if(strlen($ref) > 1)	$variant = "INDEL";	//different QC-values for SNVs and Indels		
	$as = explode(",",$alt);
	foreach($as as $a)	//check for different alleles
	{
		if($variant=="SNV" && strlen($a)>1)	$variant = "INDEL";	//different QC-values for SNVs and Indels
	}

	$o = 0;
	$r = 0;
	if($variant == "SNV" && (!is_null($idx_tumor) && !is_null($idx_normal)))
	{
		$idx_t = array_search("TU",explode(":",$row[8]));
		$idx_a = array_search("AU",explode(":",$row[8]));
		$idx_c = array_search("CU",explode(":",$row[8]));
		$idx_g = array_search("GU",explode(":",$row[8]));
		list($nuc_t,) = explode(",", explode(":", $row[$idx_tumor])[$idx_t]);
		list($nuc_a,) = explode(",", explode(":", $row[$idx_tumor])[$idx_a]);
		list($nuc_c,) = explode(",", explode(":", $row[$idx_tumor])[$idx_c]);
		list($nuc_g,) = explode(",", explode(":", $row[$idx_tumor])[$idx_g]);
		if($alt == "T") $o = $nuc_t;
		if($alt == "A") $o = $nuc_a;
		if($alt == "C") $o = $nuc_c;
		if($alt == "G") $o = $nuc_g;
		if($ref == "T") $r = $nuc_t;
		if($ref == "A") $r = $nuc_a;
		if($ref == "C") $r = $nuc_c;
		if($ref == "G") $r = $nuc_g;
	}
	else if($variant == "INDEL")
	{
		$idx_tir = array_search("TIR",explode(":",$row[8]));
		$idx_tar = array_search("TAR",explode(":",$row[8]));
		if($idx_tir!==FALSE && $idx_tar!==FALSE)
		{
			list($tir,) = explode(",", explode(":", $row[$idx_tumor])[$idx_tir]);
			list($tar,) = explode(",", explode(":", $row[$idx_tumor])[$idx_tar]);
			if(($tir+$tar)!=0)
			{
				$o = $tir;
				$r = $tar;
			}
		}
	}
	$alleles = "$r,$o";
	
	// format field: "GT  - 0/1 "column, add 
	$genotype = "1/0";
	
	$row[6] = ".";	//clear FILTER column
	$row[7] = "SOMATIC=1";	//clear INFO column
	$row[8] = "GT:AD";
	$row[$idx_tumor] = "$genotype:$alleles";
	$file_qci->addRow($row);
}

// remove normal column
if(!is_null($idx_normal))	$file_qci->removeCol($idx_normal);

// add minimal set of comments
$comments = array();
$comments[] = "#fileformat=VCF4.1";
$comments[] = "#reference=$build";
$comments[] = "#FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype. Used by QCI.'>";
$comments[] = "#FORMAT=<ID=AD,Number=.,Type=Integer,Description='Allelic depths for the ref and alt alleles in the order listed. Used by QCI.'>";
$file_qci->setComments($comments);
$vcf_qci = $parser->tempFile("_qci.vcf");
$file_qci->toTSV($vcf_qci);

//align INDELs to the left
$vcf_qci_aligned = $parser->tempFile("_qci_aligned.vcf");
$parser->exec(get_path("ngs-bits")."VcfLeftNormalize","-stream -in $vcf_qci -out $vcf_qci_aligned -ref ".genome_fasta($build), true);

//sort variants
$vcf_qci_sorted = $parser->tempFile("_qci_sorted.vcf");
//Nb: VcfStreamSort cannot be used since two vcf files (variants, indels) are concatenated and variants are not grouped by chromosome; total number of variants should be low (somatic).
$parser->exec(get_path("ngs-bits")."VcfSort","-in $vcf_qci_aligned -out $vcf_qci_sorted", true);

//zip and index output file
$parser->exec("bgzip", "-c $vcf_qci_sorted > $out", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return