<?php 
/** 
	@page an_vep
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_vep", "Variant annotation with Ensembl VEP.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addFlag("all_transcripts", "Annotate all transcripts - if unset only GENCODE basic transcripts are annotated.");

extract($parser->parse($argv));

//TODO general
//-check that custom-annotations that are now zipped are not used unzipped in other context, e.g. for CNVs
//  - OMIM needed unzipped

//VEP annotation
//TODO test plugins: dbscSNV,GeneSplicer,MaxEntScan
//TODO test: --regulatory --gene_phenotype --ccds --biotype --canonical --pubmed 
//TODO optimize speed using --buffer_size --fork
//TODO optimize VCF size by annotating only used fields
$vep_path = dirname(get_path("vep"));
$tmp = $parser->tempFile("_annotated.vcf");
$args = array();
$args[] = "-i $in"; //input
$args[] = "-o $tmp --vcf --no_stats"; //ouput
$args[] = "--species homo_sapiens --assembly {$build}"; //species
$args[] = "--offline --cache --dir_cache {$vep_path}/cache/ --fasta {$vep_path}/fasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"; //paths to data
$args[] = "--numbers --hgvs --domains"; //annotation options
$args[] = "--sift p --polyphen p"; //pathogenicity predictions
$args[] = "--af --af_gnomad --af_esp --failed 1"; //population frequencies
$args[] = "--plugin CADD,".get_path("data_folder")."/dbs/CADD/+old_1.3/whole_genome_SNVs.tsv.gz,".get_path("data_folder")."/dbs/CADD/InDels.tsv.gz"; //CADD //TODO update path
//TODO $args[] = "--plugin REVEL,".get_path("data_folder")."/dbs/REVEL/revel_all_chromosomes.csv.gz"; //REVEL
//TODO $args[] = "--plugin FATHMM_MKL,".get_path("data_folder")."/dbs/fathmm-MKL/fathmm-MKL_Current.tab.gz"; //fathmm-MKL
$args[] = "--custom ".get_path("data_folder")."/dbs/gnomAD/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF"; //genomAD genome AFs
$args[] = "--custom ".get_path("data_folder")."/dbs/RepeatMasker/RepeatMasker.bed.gz,REPEATMASKER,bed,overlap,0"; //RepeatMasker
$args[] = "--custom ".get_path("data_folder")."/dbs/ClinVar/clinvar_20180701_converted.vcf.gz,CLINVAR,vcf,exact,0,DETAILS"; //ClinVar
//TODO $args[] = "--custom ".get_path("data_folder")."/dbs/phyloP/hg19.100way.phyloP100way.bw,PHYLOP,bigwig";
if(file_exists(get_path("data_folder")."/dbs/OMIM/omim.bed.gz")) //OMIM annotation (optional because of license)
{
	$args[] = "--custom ".get_path("data_folder")."/dbs/OMIM/omim.bed.gz,OMIM,bed,overlap,0";
}
if(file_exists(get_path("data_folder")."/dbs/HGMD/HGMD_PRO_2018_2_fixed.vcf.gz")) //HGMD annotation (optional because of license)
{
	$args[] = "--custom ".get_path("data_folder")."/dbs/HGMD/HGMD_PRO_2018_2_fixed.vcf.gz,HGMD,vcf,exact,0,CLASS,MUT,GENE,PHEN";
}	
if (!$all_transcripts)
{
	$args[] = "--gencode_basic";
}
$parser->exec(get_path("vep"), implode(" ", $args), true);

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
