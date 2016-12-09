<?php
/** 
	@page db_auto_classification 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_auto_classification", "\$Rev: 868 $", "Automatic classification of high-impact variants in NGSD.");
$parser->addInt("c2_ihdb", "NGSD count cutoff (variants with a homozygous count >= cutoff are classified as class 2).", true, 30);
$parser->addFloat("c1_af", "Allele frequency cutoff (class 2 variants with AF <= cutoff are classified as class 1).", true, 0.01);
//optional
$parser->addFlag("update", "Enables update in NGSD. If unset, only a dry-run is performed!");
extract($parser->parse($argv));

//init
$db = DB::getInstance("NGSD");

//extract variants for each target region
print "#variant\t1000g\texac\tkaviar\tngsd_hom\tngsd_het\thgmd\tclinvar\tclass\tclass_new\tcomment\tNGSD-update\n";	
$db_vars = $db->executeQuery("SELECT id, chr, start, end, ref, obs, 1000g, exac, kaviar, coding FROM variant WHERE coding LIKE '%:HIGH:%'");
for ($i=0; $i<count($db_vars); ++$i)
{
	list($id, $chr, $start, $end, $ref, $obs, $tg, $exac, $kaviar, $coding) = array_values($db_vars[$i]);
	if (!contains($coding, ":HIGH:")) continue;
	
	//init
	$hom = "n/a";
	$het = "n/a";
	$hgmd = "n/a";
	$clinvar = "n/a";
	$class = $db->getValue("SELECT class FROM variant_classification WHERE variant_id='$id'", "n/a");
	$class_new = "n/a";
	$comments = array();
	
	//count hom/het samples
	$db_genos = $db->executeQuery("SELECT s.id, dv.genotype FROM detected_variant dv, processed_sample ps, sample s WHERE dv.variant_id=$id AND dv.processed_sample_id=ps.id AND ps.sample_id=s.id");
	$hom = 0;
	$het = 0;
	$samples_done = array();
	foreach($db_genos as $det)
	{
		$s_id = $det['id'];
		if (!isset($samples_done[$s_id]))
		{
			$geno = $det['genotype'];
			$het += $geno=="het";
			$hom += $geno=="hom";
			
			$samples_done[$s_id] = true;
		}
	}
	
	//classify (by NGSD)
	if ($hom>=$c2_ihdb)
	{
		$class_new = 2;
		
		//classify (by AF)
		if (is_numeric($tg) && $tg>=$c1_af) $class_new = 1;
		if (is_numeric($exac) && $exac>=$c1_af) $class_new = 1;
		if (is_numeric($kaviar) && $kaviar>=$c1_af) $class_new = 1;
	}
	
	//exclude pathogenic variants in ClinVar / HGMD variants (not those where WT is pathogenic)
	if ($class_new!="n/a")
	{
		$s = $start;
		$e = $end;
		$indel = $ref=="-" || $obs=="-" || strlen($ref)>1 || strlen($obs)>1;
		if ($indel)
		{
			$s -= 10;
			$s += 10;
		}

		//get clinvar annotation
		$clinvar = array();
		list($anno) = exec2("tabix -p vcf ".get_path("data_folder")."/dbs/ClinVar/clinvar_converted.vcf.gz $chr:$s-$e");
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
		list($anno) = exec2("tabix -p vcf ".get_path("data_folder")."/dbs/HGMD/HGMD_PRO_2016_1_fixed.vcf.gz $chr:$s-$e");
		foreach($anno as $line2)
		{
			if (contains($line2, "CLASS=DM") && contains($line2, "MUT=ALT"))
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
		
		if ($hgmd!="")
		{
			$comments[] = "SKIPPED_HGMD";
		}
		if ($clinvar!="")
		{
			$comments[] = "SKIPPED_CLINVAR";
		}
	}
	
	//comments
	if ($class!="n/a" && $class_new!="n/a")
	{
		if ($class>$class_new)
		{
			$comments[] = "DOWN";
		}
		if ($class<$class_new)
		{
			$comments[] = "SKIPPED_UP";
		}
		if ($class=="4" || $class=="5")
		{
			trigger_error("Cannot set class '$class_new' for variant '$chr:$start $ref>$obs' because it is classified as '$class'!", E_USER_ERROR);
		}
	}
	$comments = implode(", ", $comments);
	
	if (($class_new!="n/a" && $class_new!=$class) || $comments!="")
	{
		print "$chr:$start $ref>$obs\t$tg\t$exac\t$kaviar\t$hom\t$het\t$hgmd\t$clinvar\t$class\t$class_new\t$comments\t";
		
		if (!contains($comments, "SKIPPED"))
		{
			print "UPDATE";
			if($update)
			{
				$class_comment = $db->getValue("SELECT comment FROM variant_classification WHERE variant_id='$id'", "");
				$class_comment = "[$class_new] auto-classification ".date("d.m.Y")."\n\n".$class_comment;				
				$db->executeStmt("INSERT INTO variant_classification (variant_id, class, comment) VALUES ($id, '$class_new', '$class_comment') ON DUPLICATE KEY UPDATE class='$class_new',comment='$class_comment'");
			}
			else
			{
				print "(DRY-RUN)";
			}
		}
		print "\n";
	}
}

?>