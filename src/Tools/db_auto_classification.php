<?php
/** 
	@page db_auto_classification 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_auto_classification", "Automatic classification of exonic/splicing variants in NGSD.");
//optional
$parser->addFlag("non_coding", "Also classify intronic/intergenic variants.");
$parser->addFlag("update", "Enables update in NGSD. If unset, only a dry-run is performed!");
extract($parser->parse($argv));

//init
$db = DB::getInstance("NGSD");

//extract variants for each target region
$where = $non_coding ? "" : "WHERE coding LIKE '%:HIGH:%' OR coding LIKE '%:MODERATE:%' OR coding LIKE '%:LOW:%'";
$db_vars = $db->executeQuery("SELECT id, chr, start, end, ref, obs, 1000g, exac, gnomad, coding FROM variant $where ORDER BY chr ASC, start ASC");
print "##variants=".count($db_vars)."\n";
print "#number\tvariant\taf_1000g\taf_exac\taf_gnomad\taf_max\tngsd_hom\tngsd_het\tngsd_sum\thgmd\tclinvar\tclass\tclass_new\tcomment\tNGSD-update\n";	
for ($i=0; $i<count($db_vars); ++$i)
{
	list($id, $chr, $start, $end, $ref, $obs, $tg, $exac, $gnomad, $coding) = array_values($db_vars[$i]);
	
	//init
	$hom = "n/a";
	$het = "n/a";
	$hgmd = "n/a";
	$clinvar = "n/a";
	$class = $db->getValue("SELECT class FROM variant_classification WHERE variant_id='$id'", "n/a");
	$class_new = "n/a";
	$comments = array();
	
	//classify by public database frequency
	$max_af = 0.0;
	if (is_numeric($tg) && $tg>$max_af) $max_af = $tg;
	if (is_numeric($exac) && $exac>$max_af) $max_af = $exac;
	if (is_numeric($gnomad) && $gnomad>$max_af) $max_af = $gnomad;
	if ($max_af>=0.05)
	{
		$class_new = 2;
	}
	
	//classify by NGSD
	$hom = $db->getValue("SELECT count_hom FROM detected_variant_counts WHERE variant_id='$id'", "0");
	$het = $db->getValue("SELECT count_het FROM detected_variant_counts WHERE variant_id='$id'", "0");
	if ($chr!="chrMT" && $hom>=30)
	{
		$class_new = 2;
	}
	
	//get clinvar annotation
	$s = $start;
	$e = $end;
	$indel = $ref=="-" || $obs=="-" || strlen($ref)>1 || strlen($obs)>1;
	if ($indel)
	{
		$s -= 10;
		$s += 10;
	}

	$clinvar = array();
	list($anno) = exec2("tabix -p vcf ".get_path("data_folder")."/dbs/ClinVar/clinvar_20180128_converted.vcf.gz $chr:$s-$e");
	foreach($anno as $line2)
	{
		if (contains($line2, "pathogenic"))
		{
			list($c2, $s2, , $r2, $o2s, , , $sig) = explode("\t", trim($line2));
			foreach(explode(",", $o2s) as $o2)			
			{
				$indel2 = strlen($r2)>1 || strlen($o2)>1;
				
				//SNP hit
				if (!$indel && !$indel2)
				{
					$clinvar[] = "$c2:$s2 $r2>$o2s $sig";
				}
				
				//INDEL hit
				if ($indel && $indel2)
				{
					$clinvar[] = "$c2:$s2 $r2>$o2s $sig";
				}
			}
		}
	}
	$clinvar = implode(", ", $clinvar);
	
	//get HGMD annotation
	$hgmd = array();
	list($anno) = exec2("tabix -p vcf ".get_path("data_folder")."/dbs/HGMD/HGMD_PRO_2017_4_fixed.vcf.gz $chr:$s-$e");
	foreach($anno as $line2)
	{
		if (contains($line2, "CLASS=DM"))
		{
			list($c2, $s2, , $r2, $o2s, , , $sig) = explode("\t", trim($line2));
			foreach(explode(",", $o2s) as $o2)			
			{
				$indel2 = strlen($r2)>1 || strlen($o2)>1;
				
				//SNP hit
				if (!$indel && !$indel2)
				{
					$hgmd[] = "$c2:$s2 $r2>$o2s $sig";
				}
				
				//INDEL hit
				if ($indel && $indel2)
				{
					$hgmd[] = "$c2:$s2 $r2>$o2s $sig";
				}
			}
		}
	}
	$hgmd = implode(", ", $hgmd);
	
	//special handling of rare pathogenic variants according to clinvar/HGMD
	if ($class_new!="n/a" && (contains($hgmd, "CLASS=DM") || (contains($clinvar, "pathogenic") && !contains($clinvar, "conflicting"))))
	{
		$class_new = 3;
	}
	
	//comments
	if ($class!="n/a" && $class_new!="n/a")
	{
		if ($class<$class_new)
		{
			$comments[] = "SKIPPED_UP";
		}
		if ($class=="M")
		{
			$comments[] = "SKIPPED_IS_M";
		}
		if ($class=="4" || $class=="5")
		{
			$comments[] = "SKIPPED_IS_4_OR_5";
		}
	}
	$comments = implode(", ", $comments);
	
	if (($class_new!="n/a" && $class_new!=$class) || $comments!="")
	{
		print "$i\t$chr:$start $ref>$obs\t$tg\t$exac\t$gnomad\t$max_af\t$hom\t$het\t".(2*$hom+$het)."\t$hgmd\t$clinvar\t$class\t$class_new\t$comments\t";
		
		if (!contains($comments, "SKIPPED"))
		{
			print "UPDATE";
			if($update)
			{
				$class_comment = $db->getValue("SELECT comment FROM variant_classification WHERE variant_id='$id'", "");
				$class_comment = str_replace("'", "\'", $class_comment);
				$class_comment = "[$class_new] auto-classification ".date("d.m.Y")."\n\n".$class_comment;				
				$db->executeStmt("INSERT IGNORE INTO variant_classification (variant_id, class, comment) VALUES ($id, '$class_new', '$class_comment') ON DUPLICATE KEY UPDATE class='$class_new',comment='$class_comment'");
			}
			else
			{
				print " (DRY-RUN)";
			}
		}
		print "\n";
	}
}

?>