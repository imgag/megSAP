<?php 
/** 
	@page an_somatic_cancerhotspots
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_somatic_cancerhotspots", "Variant annotation with somatic cancer hotspots from cancerhotspots.org based on AA, which is why we cannot use VcfAnnotateFromVcf.");
$parser->addInfile("in",  "Input file in VCF format (variant consequence annotation with VcfAnnotateConsequence must be present).", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
extract($parser->parse($argv));

//get index of columnn in QSC header.
function index_of($cols, $name)
{
	$index = array_search($name, $cols);
	if ($index===FALSE)
	{
		trigger_error("Could not find column '$name' in QSC2 annotation. Valid column names are: ".implode(", ", array_values($cols)), E_USER_ERROR);
	}
	return $index;
}

//load cancer hotspot data
$chs_data = [];
$in_db = repository_basedir()."/data/misc/cancerhotspots/cancerhotspots_snv.tsv.gz";
$h = gzopen2($in_db, "r");
while(!gzeof($h))
{
	$line = nl_trim(gzgets($h));
	if ($line=="" || $line[0]=="#") continue;
	
	list($gene, $transcript, $pos, $ref, $alt, $total_count, $alt_count) = explode("\t", $line);
	
	$chs_data[$transcript][$ref.$pos.$alt] = implode("|",[$gene,$transcript,$pos,$ref,$alt,$total_count,$alt_count]);

}
gzclose($h);

//get annotation indices in CSQ2 field
$i_csq2_hgvsp = -1;
$i_csq2_feature = -1;
$handle_in = fopen2($in, "r");
while(!feof($handle_in))
{
	$line_in = trim(fgets($handle_in));
	if ($line_in=="") continue;
	
	//headers
	if(starts_with($line_in, "#"))
	{
		if (starts_with($line_in, "##INFO=<ID=CSQ2,"))
		{
			$cols = explode("|", substr($line_in, 0, -2));
			$i_csq2_hgvsp = index_of($cols, "HGVSp");
			$i_csq2_feature = index_of($cols, "Feature");
		}
		
		continue;
	}
	
	break; //first content line
}
fclose($handle_in);
if ($i_csq2_hgvsp==-1) trigger_error("VCF file does not contain variant consequence annotation in CSQ2 INFO entry!", E_USER_ERROR);

//perform annotation og intput file
$handle_in = fopen2($in, "r");
$handle_out = fopen2($out, "w");
while(!feof($handle_in))
{
	$line_in = trim(fgets($handle_in));
	if ($line_in=="") continue;
	
	//header
	if(starts_with($line_in, "#"))
	{
		//write INFO field for CANCERHOTSPOTS just before INFO field of CSQ annotation
		if (starts_with($line_in, "##INFO=<ID=CSQ2,"))
		{
			fwrite($handle_out,"##INFO=<ID=CANCERHOTSPOTS,Number=.,Type=String,Description=\"Cancerhotspots.org data. Format: GENE_SYMBOL|ENSEMBL_TRANSCRIPT|AA_POS|AA_REF|AA_ALT|TOTAL_COUNT|ALT_COUNT\">\n");
		}
		
		//write out header
		fwrite($handle_out, "{$line_in}\n");
		
		continue;
	}
	
	//parse variants	
	$parts = explode("\t", $line_in);
	if(count($parts)<9) trigger_error("VCF line with less than 9 fields found: {$line_in}.", E_USER_ERROR);

	//determine AA change matches for this variant
	$aa_matches = [];
	foreach(explode(";", $parts[7]) as $info_part)
	{
		if(!starts_with($info_part, "CSQ2=")) continue;
		
		foreach(explode(",", substr($info_part, 5)) as $csq_entry)
		{
			$csq_parts = explode("|", $csq_entry);
			
			//check if transcript is in cancerhotspots data
			$transcript = explode(".", $csq_parts[$i_csq2_feature].".")[0]; //remove transcript version if appended
			if (!isset($chs_data[$transcript])) continue;
			
			//skip protein changes other than single AA substitutions
			$matches = [];
			if(!preg_match("/^p.([a-zA-Z]{3})(\d+)([a-zA-Z]{3})$/", $csq_parts[$i_csq2_hgvsp], $matches)) continue;
			list(, $aa_ref, $aa_pos, $aa_alt) = $matches;
			
			//check if protein change is in cancerhotspots data
			$aa_change = aa3_to_aa1($aa_ref).$aa_pos.aa3_to_aa1($aa_alt);
			if (!isset($chs_data[$transcript][$aa_change])) continue;
			
			$aa_matches[] = $chs_data[$transcript][$aa_change];
		}
	}
	
	//add matches from cancerhotspots input file to INFO field
	if(!empty($aa_matches)) $parts[7] .= ";CANCERHOTSPOTS=" . implode("," , $aa_matches);
	
	fwrite($handle_out, implode("\t", $parts) ."\n" );
}
fclose($handle_in);
fclose($handle_out);
?>
