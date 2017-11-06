<?php


/**
	@page somatic_nearby_germline_variants
*/

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";
require_once($basedir."Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// command line arguments
$parser = new ToolBase("somatic_nearby_germline_variants", "Identify germline variants close to somatic variants.");
$parser->addInfile("somatic_var", "Somatic variants in GSvar format.", false);
$parser->addInfile("normal_bam", "Normal sample BAM file.", false);
$parser->addString("n_dna_id", "Normal sample identifier.", false);
$parser->addString("n_dna_sys", "Normal sample processing system.", false);
extract($parser->parse($argv));

$transcript_ids = array();
$transcript_type = "";
$file_som = Matrix::fromTSV($somatic_var);
$idx_cs = $file_som->getColumnIndex("coding_and_splicing");
$idx_filter = $file_som->getColumnIndex("filter");
$known_gene = get_path("data_folder")."/dbs/UCSC/knownGene.txt";
$no_transcript_id = array();

// (4a) extract transcripts and their exons for given variants
for($i=0;$i<$file_som->rows();++$i)	// get all somatic refseq-ids
{
	$row_som = $file_som->getRow($i);
//		if(!empty($filter))	continue;

	$tmp_ids = array();
	//REFSEQ transcript IDs
	//preg_match('/NM_\d+/',$row_som[$idx_cs],$tmp_ids);
	//ENSEMBL transcript IDs: ENST00000376030
	//preg_match('/NM_\d+/',$row_som[$idx_cs],$tmp_ids);
	foreach(explode(",",$row_som[$idx_cs]) as $t)
	{
		$tmp_id = explode(":",$t)[1];
		if(strpos($tmp_id,"ENST")==0)	$transcript_type = "Ensembl";
		else if(strpos($tmp_id,"NM_")==0)	$transcript_type = "RefSeq";
		else if(!empty($tmp_id))	trigger_error("Unkown transcript information $tmp_id.",E_USER_ERROR);
		$tmp_ids[] = $tmp_id;
	}

	if(count($tmp_ids) == 0)
	{
		$no_transcript_id[] = $i;
		continue;
	}

	foreach($tmp_ids as $tmp_id)
	{
		if(!isset($transcript_ids[$tmp_id]))
		{
			$transcript_ids[$tmp_id] = array();
			$transcript_ids[$tmp_id]['variant_rows'] = array();
		}
		$transcript_ids[$tmp_id]['variant_rows'][] = $i;
	}
}
if(!empty($no_transcript_id))	trigger_error ("No transcript ID was identified in rows: ".implode(", ", $no_transcript_id).". Variants were skipped.",E_USER_NOTICE);

//convert transcript-ids to UCSC-ids
$transcript_file = "";
if($transcript_type=="Ensembl")	$transcript_file = get_path("data_folder")."/dbs/UCSC/knownToEnsembl.tsv";
if($transcript_type=="RefSeq")	$transcript_file = get_path("data_folder")."/dbs/UCSC/knownToRefSeq.tsv";

//$kgxref_file = get_path("data_folder")."/dbs/UCSC/kgXref.txt";
$ucsc2transcript = array();

$handle = fopen($transcript_file, "r");
if($handle)
{
	while(($buffer=fgets($handle)) !== FALSE)
	{
		$row = explode("\t", trim($buffer));

		if(array_key_exists($row[1], $transcript_ids))	$ucsc2transcript[$row[0]] = $row[1];
	}
	fclose($handle);
}
else	trigger_error("Could not find file $transcript_file.",E_USER_ERROR);

// extract cDNA information from UCSC-ids
$handle = fopen($known_gene, "r");
if($handle)
{
	while(($buffer=fgets($handle)) !== FALSE)
	{
		$row = explode("\t", trim($buffer));

		if(array_key_exists($row[0], $ucsc2transcript))
		{
			$exon_start = explode(",", $row[8]);
			$exon_end = explode(",", $row[9]);
			$transcript_ids[$ucsc2transcript[$row[0]]]['strand'] = $row[2];
			$transcript_ids[$ucsc2transcript[$row[0]]]['cdsStart'] = $row[5];
			$transcript_ids[$ucsc2transcript[$row[0]]]['cdsEnd'] = $row[6];
			$transcript_ids[$ucsc2transcript[$row[0]]]['exons'] = array();
			for($i=0;$i<count($exon_start);++$i)
			{
				if(empty($exon_start[$i]))	continue;
				$transcript_ids[$ucsc2transcript[$row[0]]]['exons'][] = array('start'=>$exon_start[$i], 'end'=>$exon_end[$i]);
			}
		}
	}
	fclose($handle);
}
else trigger_error("Could not open file $known_gene.",E_USER_ERROR);

// (4b) generate region filter to filter germline variant list
$max_distance = 40;
$regions = new Matrix();
$splicing = array();
$no_exons = array();
foreach($transcript_ids as $id => $transcript_id)
{
	if(!isset($transcript_id['exons']))
	{
		$no_exons[] = $id;
		continue;
	}
	foreach($transcript_id['variant_rows'] as $variant_row)
	{
		$somatic_variant = $file_som->getRow($variant_row);
		$region = array('',-1,-1,-1,-1);
		$region[0] = $somatic_variant[0];
		// get current exon
		$current_exon = -1;
		foreach($transcript_id['exons'] as $exon => $r)
		{
			// exonic variant
			if($somatic_variant[1]>=$r['start'] && $somatic_variant[2]<=$r['end'])	$current_exon = $exon;
		}
		if($current_exon==-1)
		{
			// splicing variant, uncomment following line if region of splicing variants is relevant, TODO: distinguish between splicing, intronic, intergenic, etc.
			// $regions->addRow(array($somatic_variant[0],$somatic_variant[1]-1,$somatic_variant[2],"",$id));
			$splicing[] = $variant_row;
			continue;
		}
		// get start
		$found = false;
		$tmp_start = $somatic_variant[1];
		$tmp_diff = $max_distance;
		$tmp_exon = $current_exon;
		while($found==FALSE)
		{
			$tmp_diff -= ($tmp_start - $transcript_id['exons'][$tmp_exon]['start']);
			if($tmp_diff>0)
			{
				--$tmp_exon;
				if(!isset($transcript_id['exons'][$tmp_exon]['start']))
				{
					trigger_error("Could not determine start point for variant in row $variant_row with starting position ".$somatic_variant[1].", using ".$transcript_id['exons'][$tmp_exon+1]['start']." as start point (first exon).",E_USER_NOTICE);
					$region[1] = $transcript_id['exons'][$tmp_exon+1]['start'];
					$region[3] = $tmp_exon+1;
					$found=TRUE;
					continue;
				}
				$tmp_start = $transcript_id['exons'][$tmp_exon]['end'];
			}
			else if($tmp_diff<=0)
			{
				$region[1] = $transcript_id['exons'][$tmp_exon]['start'] - $tmp_diff - 1;
				$region[3] = $tmp_exon;
				$found=TRUE;
			}
		}
		// get end
		$found = false;
		$tmp_end = $somatic_variant[2];
		$tmp_diff = $max_distance;
		$tmp_exon = $current_exon;
		while($found==FALSE)
		{
			$tmp_diff -= ($transcript_id['exons'][$tmp_exon]['end'] - $tmp_end);
			if($tmp_diff>0)
			{
				++$tmp_exon;
				if(!isset($transcript_id['exons'][$tmp_exon]['start']))
				{
					trigger_error("Could not determine endpoint for variant in row $variant_row with starting position ".$somatic_variant[1].", using last exon as end point.",E_USER_NOTICE);
					$region[2] = $transcript_id['exons'][$tmp_exon-1]['end'];
					$region[4] = $tmp_exon-1;
					$found=TRUE;
					continue;
				}
				$tmp_end = $transcript_id['exons'][$tmp_exon]['start'];
			}
			else if($tmp_diff<=0)
			{
				$region[2] = $transcript_id['exons'][$tmp_exon]['end'] + $tmp_diff;
				$region[4] = $tmp_exon;
				$found=TRUE;
			}
		}
		// extract all exons from start to end
		$variant_id = implode("_", array($somatic_variant[0],$somatic_variant[1],$somatic_variant[2],$somatic_variant[3],$somatic_variant[4]));
		for($exon=$region[3];$exon<=$region[4];++$exon)
		{
			$tmp_start = $transcript_id['exons'][$exon]['start'];
			$tmp_end = $transcript_id['exons'][$exon]['end'];
			if($tmp_start < $region[1])	$tmp_start = $region[1];
			if($tmp_end > $region[2])	$tmp_end = $region[2];
			$regions->addRow(array($region[0],$tmp_start, $tmp_end,$exon+1,$id,$variant_id));
		}
	}
}
if(!empty($no_exons))	trigger_error("No exons availabe for ".count($no_exons)." IDs: ".implode(", ",$no_exons),E_USER_NOTICE);
if(!empty($splicing))	trigger_error("Skipped splicing variants in the following rows: ".implode(", ",array_unique($splicing)).". Variants were skipped.", E_USER_NOTICE);
$out_reg = $parser->tempFile("_reg.bed");
$regions->toTSV($out_reg);
// $regions->toTSV( $o_folder.basename($nor_bam, ".bam").".reg");

// (4c) find reference variants around coding variants and generate filtered reference variant file
$tmp_folder = $parser->tempFolder("analyze");
$nor_vcffile = $tmp_folder."/".$n_dna_id.".vcf.gz";
$nor_varfile = $tmp_folder."/".$n_dna_id.".GSvar";
$extras_vc = "-target $out_reg ";
$n_dna_system = load_system($n_dna_sys);
$build = $n_dna_system['build'];
$parser->execTool("NGS/vc_freebayes.php", "-bam $normal_bam -out $nor_vcffile -build ". $build." $extras_vc");
$parser->execTool("Pipelines/annotate.php", "-out_name $n_dna_id -out_folder $tmp_folder -system $n_dna_sys -vcf $nor_vcffile");

// add new col that contains information about surrounding variants
$new_col = array();
$ref = Matrix::fromTSV($nor_varfile);
for($i=0;$i<$file_som->rows();++$i)	// get all somatic refseq-ids
{
	$row_som = $file_som->getRow($i);
	$variant_id = implode("_", array($row_som[0],$row_som[1],$row_som[2],$row_som[3],$row_som[4]));

	//identify regions that belong to this variant
	$target_regions = array();
	for($j=0;$j<$regions->rows();++$j)
	{
		$row_reg = $regions->getRow($j);
		if($row_reg[5]==$variant_id)
		{
			$target_regions[] = array($row_reg[0],$row_reg[1]+1,$row_reg[2]);
		}
	}

	// check if germline variants are in these regions
	$has_germline_variant = false;
	for($j=0;$j<$ref->rows();++$j)
	{
		$row_ref = $ref->getRow($j);

		foreach($target_regions as $target_region)
		{
			if($row_ref[0]==$target_region[0])
			{
				if(range_overlap($row_ref[1],$row_ref[2],$target_region[1],$target_region[2]))	$has_germline_variant = true;
			}
		}
	}

	if($has_germline_variant)	$new_col[] = "yes";
	else	$new_col[] = "no";
}
$file_som->addCol($new_col, "nearby_germline_variant","There are germline variants surrounding this variant (+/- $max_distance bp, new feature - use with caution, haplotype not considered).");
$file_som->toTSV($somatic_var);