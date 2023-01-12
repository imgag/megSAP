<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//format AF
function format_af($an, $af)
{
	//gnomAD 3.1.2 contains 75000 genomes, i.e. 150000 alleles. If the allele number is below 1000, we return 0 because the region is not reliable
	if ($an<1000) return "0.0";
	
	//round AF to 5 digits
	$af_before = $af;
	$af = number_format($af, 5);
	
	//make sure that the rounding does not set AF=0 although it is greater than 0
	if ($af==0.0 && $af_before>0) $af="0.00001";
	
	return $af;
}

$handle = fopen2("php://stdin", "r");
while($line = fgets($handle))
{
	$line = trim($line);
	if ($line=='') continue;

	if($line[0]=='#')
	{
		if (in_array('-header', $argv))
		{
			if($line[1]=='#')
			{
				//comment section

				//skip INFO ID headers (not present in output file)
				if (starts_with($line, '##INFO=<ID=')) continue;

				//print comment line
				print $line."\n";
			}
			else
			{

				//actual header line
				// add INFO column header:
				print "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in samples\">\n";
				print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency in samples\">\n";
				print "##INFO=<ID=AFR_AF,Number=A,Type=Float,Description=\"Alternate allele frequency in AFR samples\">\n";
				print "##INFO=<ID=AMR_AF,Number=A,Type=Float,Description=\"Alternate allele frequency in AMR samples\">\n";
				print "##INFO=<ID=EAS_AF,Number=A,Type=Float,Description=\"Alternate allele frequency in EAS samples\">\n";
				print "##INFO=<ID=NFE_AF,Number=A,Type=Float,Description=\"Alternate allele frequency in NFE samples\">\n";
				print "##INFO=<ID=SAS_AF,Number=A,Type=Float,Description=\"Alternate allele frequency in SAS samples\">\n";
				print "##INFO=<ID=AC,Number=A,Type=Float,Description=\"Alternate allele count\">\n";
				print "##INFO=<ID=Hom,Number=A,Type=Integer,Description=\"Count of homozygous individuals in samples\">\n";
				print "##INFO=<ID=Hemi,Number=A,Type=Integer,Description=\"Alternate allele count for male samples\">\n";
				print "##INFO=<ID=Het,Number=A,Type=Integer,Description=\"Count of heterozygous alleles in all samples.\">\n";
				print "##INFO=<ID=Wt,Number=A,Type=Integer,Description=\"Count of wildtype alleles in all samples.\">\n";
				
				//print header line
				print $line."\n";
			}
		}
		continue;
	}
	
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = explode("\t", $line);
		
	//parse info entries
	$nonpar = false;
	$ac = null;
	$af = null;
	$an = null;
	$hom = null;
	$ac_xy = null;
	$af_afr = null;
	$af_amr = null;
	$af_eas = null;
	$af_nfe = null;
	$af_sas = null;
	foreach(explode(';', $info) as $entry)
	{
		$sep_idx = strpos($entry, '=');
		if ($sep_idx===false)
		{
			if ($entry=='nonpar') $nonpar = true;
		}
		else
		{
			$key = substr($entry, 0, $sep_idx);
			if ($key=='AC')
			{
				$ac = substr($entry, $sep_idx+1);
			}
			else if ($key=='AF')
			{
				$af = substr($entry, $sep_idx+1);
			}
			else if ($key=='AN')
			{
				$an = substr($entry, $sep_idx+1);
			}
			else if ($key=='nhomalt')
			{
				$hom = substr($entry, $sep_idx+1);
			}
			else if ($key=='AC_XY')
			{
				$ac_xy = substr($entry, $sep_idx+1);
			}
			else if ($key=='AF_afr')
			{
				$af_afr = substr($entry, $sep_idx+1);
			}
			else if ($key=='AF_amr')
			{
				$af_amr = substr($entry, $sep_idx+1);
			}
			else if ($key=='AF_eas')
			{
				$af_eas = substr($entry, $sep_idx+1);
			}
			else if ($key=='AF_nfe')
			{
				$af_nfe = substr($entry, $sep_idx+1);
			}
			else if ($key=='AF_sas')
			{
				$af_sas = substr($entry, $sep_idx+1);
			}
		}
	}
	
	$info_new = array();
	$info_new[] = 'AN='.$an;
	if ($chr=='chrX' && $nonpar && ! is_null($ac_xy))
	{
		$info_new[] = 'Hemi='.$ac_xy;
	}
	else if ($chr=='chrY')
	{
		$info_new[] = 'Hemi='.$ac; // full AC for chrY variants as they are always hemizygote
	}
	else
	{
		$info_new[] = 'Hemi=.';
	}
	
	$info_new[] = 'Hom='.$hom;
	$info_new[] = 'Het='.($ac - 2*$hom);
	$info_new[] = 'Wt='.($an-$ac);
	$info_new[] = 'AC='.$ac;
	$info_new[] = 'AF='.format_af($an, $af);
	$info_new[] = 'AFR_AF='.format_af($an, $af_afr);
	$info_new[] = 'AMR_AF='.format_af($an, $af_amr);
	$info_new[] = 'EAS_AF='.format_af($an, $af_eas);
	$info_new[] = 'NFE_AF='.format_af($an, $af_nfe);
	$info_new[] = 'SAS_AF='.format_af($an, $af_sas);

	print "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t".implode(';', $info_new)."\n";
}
fclose($handle);

?>