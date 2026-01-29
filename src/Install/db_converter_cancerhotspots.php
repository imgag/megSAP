<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//convert gene to approved symbol
function gene2approved($gene)
{
	static $approved = [];
	
	$gene = trim($gene);
	
	if (!isset($approved[$gene]))
	{
		list($stdout) = exec2("echo \"{$gene}\" | GenesToApproved | cut -f1");
		$approved[$gene] = trim(implode("", $stdout));
	}

	return $approved[$gene];
}

//convert gene to best transcript
function gene2transcript($gene)
{
	global $gene2trans_mane;
	global $gene2trans_best;
	
	static $best = [];
	
	$gene = trim($gene);
	
	if (!isset($best[$gene]))
	{
		$trans = "";
		
		//check if there is a MANE select transcript
		list($stdout) = exec2("grep \"{$gene}\" {$gene2trans_mane}", false);
		foreach($stdout as $line)
		{
			$line = trim($line);
			if (!starts_with($line, $gene."\t")) continue;
			$trans = trim(explode("\t", $line)[1]);
		}
		
		//if no MANE select, just take the best transcript
		if ($trans=="")
		{
			list($stdout) = exec2("grep \"{$gene}\" {$gene2trans_best}", false);
			foreach($stdout as $line)
			{
				$line = trim($line);
				if (!starts_with($line, $gene."\t")) continue;
				$trans = trim(explode("\t", $line)[1]);
			}
		}
		
		$best[$gene] = $trans;
	}

	return $best[$gene];
}

// parse command line arguments
$parser = new ToolBase("db_converter_cancerhotspots", "Converts the cancer hotspots MAF file to TSV.");
$parser->addInfile("in", "Input file in MAF.GZ format.", false);
$parser->addOutfile("out", "Output file in TSV format.", false);
extract($parser->parse($argv));

//parse input (extract variants and other relevant information)
print "parsing variants from MAF file\n";
$c_all = 0;
$c_skip_not_snp = 0;
$c_skip_not_hg19 = 0;
$var_data = [];
$genes = [];
$handle = gzopen2($in, "r");
while(!feof($handle))
{
	$line = trim(fgets($handle));
	if($line=="") continue;

	//header > check version
	if($line[0]=="#")
	{
		if ($line!="#version 2.4") trigger_error("Ivalid version '{$line}'! 2.4 expected!", E_USER_ERROR);
		continue;
	}
	
	//variant entries
	$parts = explode("\t", $line);
	
	//gene
	$gene = $parts[0];
	$gene = gene2approved($gene);
	$genes[$gene] = true;
	
	++$c_all;
	
	//check variant is in GRCh37 genome
	if ($parts[3]!="GRCh37")
	{
		++$c_skip_not_hg19;
		continue;
	}
	
	//check variant type
	$type = $parts[9];
	if ($type!="SNP")
	{
		++$c_skip_not_snp;
		continue;
	}
	
	//variant
	$chr = "chr".$parts[4];
	$pos = $parts[5];
	$ref = $parts[10];
	$alt1 = $parts[11];
	$alt2 = $parts[12];
	if ($alt1!=$ref && $alt1!=$alt2) print "  UNEXPECTED: {$chr}:{$pos} {$ref}>{$alt1}/{$alt2}\n";
	
	//AA data
	$aa_pos = $parts[136];
	$aa_ref = $parts[138];
	$aa_alt = $parts[139];
	
	//add variant if not already contained
	$tag = "{$chr}\t{$pos}\t.\t{$ref}\t{$alt2}";
	if (!isset($var_data[$tag])) $var_data[$tag] = [0, "SOURCE={$gene}_{$chr}_{$pos}_{$ref}_{$alt2}_{$aa_ref}_{$aa_pos}_{$aa_alt}"];
	$var_data[$tag][0] += 1;
}
fclose($handle);
print "  variants skipped because not SNP: {$c_skip_not_snp}/{$c_all}\n";
print "  variants skipped because reference genome is not HG19: {$c_skip_not_hg19}/{$c_all}\n";

//prepare gene>transcript tables
$gene2trans_in = $parser->tempFile(".tsv");
file_put_contents($gene2trans_in, implode("\n", array_keys($genes)));
$gene2trans_mane = $parser->tempFile(".tsv");
exec2("GenesToTranscripts -in {$gene2trans_in} -out {$gene2trans_mane} -mode mane_select");
$gene2trans_best = $parser->tempFile(".tsv");
exec2("GenesToTranscripts -in {$gene2trans_in} -out {$gene2trans_best} -mode best");

//write VCF
$vcf = [];
$vcf[] = "##fileformat=VCFv4.3";
$vcf[] = "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Information from cancerhotspots MAF file.\">";
$vcf[] = "##INFO=<ID=COUNT,Number=1,Type=Integer,Description=\"Variant count in cancerhotspots MAF file.\">";
$vcf[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
foreach($var_data as $tag => $data)
{
	$vcf[] = "{$tag}\t30\tPASS\tCOUNT=".$data[0].";".$data[1];
}
$vcf_tmp = $parser->tempFile(".vcf");
file_put_contents($vcf_tmp, implode("\n", $vcf));

//lift-over VCF
print "lifting variants from HG19 to HG38\n";
$vcf_tmp2 = $parser->tempFile(".vcf");
list ($stdout, $stderr) = exec2("/home/ukt.ad.local/ahsturm1/.local/share/pipx/venvs/crossmap/bin/CrossMap vcf /mnt/storage1/share/opt/liftOver/hg19ToHg38.over.chain.gz {$vcf_tmp} /tmp/local_ngs_data_GRCh38/GRCh38.fa {$vcf_tmp2}"); //absolute paths are or since this script is only executed once at IMGAG
foreach($stderr as $line)
{
	$line = trim($line);
	if (contains($line, "Total entries:"))
	{
		print "  total variant entries: ".trim(explode("Total entries:", $line)[1])."\n"; 
	}
	if (contains($line, "Failed to map:"))
	{
		print "  failed to lift: ".trim(explode("Failed to map:", $line)[1])."\n"; 
	}
}

//annotate consequence
print "annotating variant consequences\n";
$vcf_tmp3 = $parser->tempFile(".vcf");
list ($stdout, $stderr) = exec2("VcfAnnotateConsequence -in {$vcf_tmp2} -gff /mnt/storage2/megSAP/data/dbs/Ensembl/Homo_sapiens.GRCh38.115.gff3 -out {$vcf_tmp3}"); //absolute paths are ok since this script is only executed once at IMGAG
foreach($stderr as $line)
{
	$line = trim($line);
	if (contains($line, "invalid variants."))
	{
		print "  ".strtolower($line)."\n"; 
	}
}

//count AA changes
print "processing AA changes\n";
$aa_data = [];
$c_no_csq = 0;
$c_no_transcript = 0;
$c_no_aa_change_match = 0;
foreach(file($vcf_tmp3) as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	if (!contains($line, "CSQ="))
	{
		++$c_no_csq;
		continue;
	}
	
	list($chr, $pos, , $ref, $alt, , , $info) = explode("\t", $line);
	
	//parse INFO entries
	$source = "";
	$gene = "";
	$count = "";
	$csq = "";
	foreach(explode(";", $info) as $entry)
	{
		if (starts_with($entry, "SOURCE="))
		{
			$source = explode("=", $entry, 2)[1];
			list($gene) = explode("_", $source, 2);
		}
		if (starts_with($entry, "COUNT="))
		{
			$count = explode("=", $entry, 2)[1];
		}
		if (starts_with($entry, "CSQ="))
		{
			$csq = explode("=", $entry, 2)[1];
		}
	}
	
	//get transcript
	$transcript = gene2transcript($gene);
	if ($transcript=="")
	{
		++$c_no_transcript;
		continue;
	}
	
	//init
	if (!isset($aa_data[$gene]))
	{
		$aa_data[$gene] = [$transcript, [], []];
	}
	
	//count
	foreach(explode(",", $csq) as $csq_part) //part example: A|missense_variant|MODERATE|BRAF|HGNC:1097|ENST00000496384.7|Transcript|13/19||c.1523C>T|p.Thr508Ile
	{
		if (!contains($csq_part, $transcript)) continue;
		
		$parts = explode("|", $csq_part);
		
		//skip protein changes other than single AA substitutions
		$matches = [];
		if(!preg_match("/^p.([a-zA-Z]{3})(\d+)([a-zA-Z]{3})$/", $parts[10], $matches))
		{
			++$c_no_aa_change_match ;
			continue;
		}
		list(, $aa_ref, $aa_pos, $aa_alt) = $matches;
		$aa_ref = aa3_to_aa1($aa_ref);
		$aa_alt = aa3_to_aa1($aa_alt);
		
		@$aa_data[$gene][1][$aa_ref.$aa_pos] += $count;
		@$aa_data[$gene][2][$aa_ref.$aa_pos.$aa_alt] += $count;
		print "$chr $pos $ref $alt // ".$parts[10]."\n";
	}
}
print "  variants skipped because not consequene annotation was available: {$c_no_csq}\n";
print "  variants skipped because gene could not be converted to transcript: {$c_no_transcript}\n";
print "  variants skipped because AA change could not be parsed: {$c_no_aa_change_match}\n";

//write output
print "writing output\n";
$out_fp = fopen2($out, "w");
fwrite($out_fp, "#GENE\tENSEMBL_TRANSCRIPT_ID\tAA_POS\tAA_REF\tAA_ALT\tCOUNT_TOTAL\tCOUNT_ALT\n");
foreach($aa_data as $gene => list($transcript, $counts_refpos, $counts_refposalt))
{
	foreach($counts_refposalt as $refposalt => $count)
	{
		$ref = $refposalt[0];
		$pos = substr($refposalt, 1, -1);
		$alt = $refposalt[strlen($refposalt)-1];
		
		$count_total = $counts_refpos[$ref.$pos];
		if ($count_total<5) continue;
		
		fwrite($out_fp, implode("\t", [$gene, $transcript, $pos, $ref, $alt, $count_total, $count])."\n");
	}
}

?>