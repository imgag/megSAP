<?php 
/** 
	@page an_snpeff
	
	@todo Look into the SSR annotation (Variant Suspect Reason Codes)
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_snpeff", "Variant annotation with SnpEff.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
//optional
$parser->addInt("thres", "Splicing region size used for annotation (flanking the exons).", true, 20);
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFlag("no_updown", "Do not annotate upsteam/downstream changes.");
extract($parser->parse($argv));

if ($build=="GRCh37")
{
	$build = "GRCh37.75";
}

//dbSNP annotation
$pipeline = array();
$pipeline[] =  array(get_path("SnpSift"), "annotate -tabix -noLog -noInfo -id ".get_path("data_folder")."/dbs/dbSNP/dbsnp_b150.vcf.gz $in");

//workaround for crash in vcfannotate when input is an empty variant list
$empty = true;
$handle = fopen($in, "r");
while(!feof($handle))
{
	$line = trim(fgets($handle));
	if ($line!="" && $line[0]!="#")
	{
		$empty = false;
		break;
	}
}
if (!$empty)
{
	//OMIM annotation (optional because of license)
	$db_file = get_path("data_folder")."/dbs/OMIM/omim.bed";
	if(file_exists($db_file))
	{
		$pipeline[] =  array(get_path("ngs-bits")."VcfAnnotateFromBed", "-bed $db_file -name OMIM");
	}

	//RepeatMasker annotation
	$pipeline[] =  array(get_path("ngs-bits")."VcfAnnotateFromBed", "-bed ".get_path("data_folder")."/dbs/RepeatMasker/RepeatMasker.bed -name REPEATMASKER");
}

//SnpEff annotation
$args = array();
if ($no_updown)
{
	$args[] = "-no-downstream ";
	$args[] = "-no-upstream ";
}
$pipeline[] =  array(get_path("SnpEff"), "eff -noLog -noStats -noInteraction -spliceRegionIntronMax $thres $build ".implode(" ",$args));

//dbNSFP annotation
$cols = array("phyloP100way_vertebrate","SIFT_pred","MetaLR_pred","FATHMM_pred","Polyphen2_HVAR_pred","CADD_phred");
$pipeline[] =  array(get_path("SnpSift"), "dbnsfp -noLog -db ".get_path("data_folder")."/dbs/dbNSFP/dbNSFPv2.9.3.txt.gz -f ".implode(",",$cols)." -");

//HGMD annotation (optional because of license)
$db_file = get_path("data_folder")."/dbs/HGMD/HGMD_PRO_2017_4_fixed.vcf";
if(file_exists($db_file))
{
	$pipeline[] =  array(get_path("SnpSift"), "annotate -mem -sorted -noLog -noId -name HGMD_ -info ID,CLASS,MUT,GENE,PHEN $db_file");
}

//gnomAD annotation
$pipeline[] =  array(get_path("SnpSift"), "annotate -tabix -noLog -noId -name GNOMAD_ -info AF ".get_path("data_folder")."/dbs/gnomAD/gnomAD_r2.0.2.vcf.gz");
$pipeline[] =  array(get_path("SnpSift"), "annotate -tabix -noLog -noId -name GNOMAD_GENOME_ -info AF ".get_path("data_folder")."/dbs/gnomAD/gnomAD_genome_r2.0.2.vcf.gz");

//1000g annotation
$pipeline[] =  array(get_path("SnpSift"), "annotate -tabix -noLog -noId -name T1000GP_ -info AF ".get_path("data_folder")."/dbs/1000G/1000g_v5b.vcf.gz");

//ExAC annotation
$pipeline[] =  array(get_path("SnpSift"), "annotate -tabix -noLog -noId -name EXAC_ -info AF,AC_Hom,Hom_NFE,Hom_AFR,AF_AFR,AF_AMR,AF_EAS,AF_NFE,AF_SAS ".get_path("data_folder")."/dbs/ExAC/ExAC_r0.3.1.vcf.gz");

//ClinVar annotation
$pipeline[] =  array(get_path("SnpSift"), "annotate -mem -sorted -noLog -noId -name CLINVAR_ -info SIG,ACC,DISEASE ".get_path("data_folder")."/dbs/ClinVar/clinvar_20180128_converted.vcf");

//COSMIC annotation (optional because of license)
$db_file = get_path("data_folder")."/dbs/COSMIC/cosmic.vcf.gz";
if(file_exists($db_file))
{
	$pipeline[] =  array(get_path("SnpSift"), "annotate -tabix -noLog -noId -name COSMIC_ -info ID $db_file");
}

//execute pipeline
$tmp = $parser->tempFile("_annotated.vcf");
$pipeline[count($pipeline)-1][1] .= " > $tmp";
$parser->execPipeline($pipeline, "annotation");

//broken info field headers (otherwise check_vcf fails)
$invalid_num_headers = array("RO","GTI","NS","SRF","NUMALT","DP","QR","SRR","SRP","PRO","EPPR","DPB","PQR","RPPR","MQMR","ODDS","AN","PAIREDR", //FreeBayes and VcfLib
							 "HGMD_GENE", "HGMD_CLASS", "HGMD_MUT", "HGMD_PHEN", //HGMD
							 );
//missing info field headers (otherwise check_vcf fails)
$comments = array();
$comments[] = "##INFO=<ID=HGMD_ID,Number=.,Type=String,Description=\"HGMD identifier(s)\">\n";
$comments[] = "##INFO=<ID=EXAC_AF_AFR,Number=1,Type=Float,Description=\"ExAC AFR subpopulation allele frequency.\">\n";
$comments[] = "##INFO=<ID=EXAC_AF_AMR,Number=1,Type=Float,Description=\"ExAC AMR subpopulation allele frequency.\">\n";
$comments[] = "##INFO=<ID=EXAC_AF_EAS,Number=1,Type=Float,Description=\"ExAC EAS subpopulation allele frequency.\">\n";
$comments[] = "##INFO=<ID=EXAC_AF_NFE,Number=1,Type=Float,Description=\"ExAC NFE subpopulation allele frequency.\">\n";
$comments[] = "##INFO=<ID=EXAC_AF_SAS,Number=1,Type=Float,Description=\"ExAC SAS subpopulation allele frequency.\">\n";


$handle1 = fopen($tmp, "r");
$handle2 = fopen($out, "w");
$comments_written = false;
while(!feof($handle1))
{
	$line = fgets($handle1);
	
	if (starts_with($line, "##")) //handle comments
	{
		if ($comments_written)
		{
			trigger_error("Invalid VCF header - all comment lines must be at the beginning. This line is not:\n$line", E_USER_ERROR);
		}
		
		//fix broken info field headers
		foreach($invalid_num_headers as $header)
		{
			if (starts_with($line, "##INFO=<ID={$header},"))
			{
				$line = str_replace(",Number=1,", ",Number=.,", $line);
			}
		}
		
		$comments[] = $line;
	}
	else if (starts_with($line, "#")) //handle header line
	{
		//sort and write comments headers before header line
		$comments = sort_vcf_comments($comments);
		foreach($comments as $comment)
		{
			fwrite($handle2, $comment);
		}
		$comments_written = true;
		
		fwrite($handle2, $line);
	}
	else //handle content lines
	{
		fwrite($handle2, $line);
	}
}
fclose($handle1);
fclose($handle2);

?>
