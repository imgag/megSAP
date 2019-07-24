<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$ngsbits = get_path("ngs-bits");

function gene2coords($gene)
{
	global $ngsbits;
	
	list($stdout, $stderr) = exec2("echo '{$gene}' | $ngsbits/GenesToBed -source ensembl -mode gene");
	
	//error
	$error = trim(implode("", $stderr));
	if ($error!="")
	{
		return [null, null, null, null, "Could not convert gene '{$gene}' to coordinates: {$error}"];
	}
	
	//convert
	$start = null;
	$end = null;
	foreach($stdout as $line2)
	{
		$line2 = trim($line2);
		if ($line2=="") continue;
		list($chr, $s, $e, $gene_approved) = explode("\t", $line2);
		if (is_null($start))
		{
			$start = $s;
			$end = $e;
		}
		else
		{
			$start = min($s, $start);
			$end = max($e, $end);
		}
	}
	
	return [$chr, $start, $end, $gene_approved, null];
}

$in = fopen2("php://stdin", "r");
while(!feof($in))
{
	$line = trim(fgets($in));
	
	/* DELETION
		`disease` varchar(125) COMMENT 'HGMD phenotype as reported in cited reference'
		`gene` varchar(15) COMMENT 'HGMD gene symbol (from MARKNAME table)'
		`descr` varchar(75) COMMENT 'Narrative description of variant'
		`cdna` char(1) COMMENT 'Reported at genomic (g) or cDNA (c) level'
		`tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COMMENT 'HGMD variant class'
		`author` varchar(50) COMMENT 'First author of cited reference'
		`journal` varchar(10) COMMENT 'HGMD journal or LSDB abbreviation (from JOURNAL table)'
		`fullname` varchar(55) COMMENT 'NLM journal title abbreviation (from JOURNAL table)'
		`vol` varchar(15) COMMENT 'Volume of cited reference'
		`page` varchar(25) COMMENT 'First page of cited reference'
		`year` year(4) COMMENT 'Publication year of cited reference'
		`pmid` varchar(8) COMMENT 'Pubmed ID of cited reference'
		`comments` text COMMENT 'HGMD curator comments'
		`acc_num` varchar(10) NOT NULL DEFAULT '' COMMENT 'HGMD mutation accession number'
		`new_date` date COMMENT 'Date entered into HGMD'
	*/
	if (starts_with($line, "INSERT INTO `grosdel` VALUES"))
	{
		$line = substr($line, 31, -3);
		$entries = explode("'),('", $line);
		foreach($entries as $entry)
		{
			$entry = "'".$entry."'";
			//print $entry."\n";
			$parts = str_getcsv($entry, ',', '\'');
			//print_r($parts);
			if (count($parts)!=15) trigger_error("Expected 15 parts, but got ".count($parts)."!", E_USER_ERROR);
			
			list($disease, $gene, $desc, $cdna, $tag, $author, $journal, $fullname, $vol, $page, $year, $pmid, $comments, $acc_num, $new_date) = $parts;
			
			$id = trim($acc_num);
			$class = trim($tag);
			$disease = trim($disease);
			$disease = strtr($disease, '"', '\'');
			$desc = trim($desc);
			$desc = strtr($desc, '"', '\'');
			$pmid = trim($pmid);
			if ($pmid=="NULL") $pmid = "";
			
			//convert gene to coordinates
			list($chr, $start, $end, $gene_approved, $error) = gene2coords($gene);
			if (!is_null($error))
			{
				trigger_error($error." ($id)", E_USER_WARNING);
				continue;
			}
			print "{$chr}\t{$start}\t{$end}\t{$id} [GENE={$gene_approved} CLASS={$class} TYPE=deletion DISEASE=\"{$disease}\" DESC=\"{$desc}\" PMID={$pmid}]\n";
		}
	}

	/* INSERTION/DUPLICATION
		 `disease` varchar(125) COMMENT 'HGMD phenotype as reported in cited reference'
		 `gene` varchar(15) COMMENT 'HGMD gene symbol (from MARKNAME table)'
		 `type` char(1) COMMENT 'Insertion (I) or duplication (D)'
		 `descr` varchar(75) COMMENT 'Narrative description of variant'
		 `cdna` char(1) COMMENT 'Reported at genomic (g) or cDNA (c) level'
		 `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COMMENT 'HGMD variant class'
		 `author` varchar(50) COMMENT 'First author of cited reference'
		 `journal` varchar(10) COMMENT 'HGMD journal or LSDB abbreviation (from JOURNAL table)'
		 `fullname` varchar(55) COMMENT 'NLM journal title abbreviation (from JOURNAL table)'
		 `vol` varchar(15) COMMENT 'Volume of cited reference'
		 `page` varchar(25) COMMENT 'First page of cited reference'
		 `year` year(4) COMMENT 'Publication year of cited reference'
		 `pmid` varchar(8) COMMENT 'Pubmed ID of cited reference'
		 `comments` text COMMENT 'HGMD curator comments'
		 `acc_num` varchar(10) NOT NULL DEFAULT '' COMMENT 'HGMD mutation accession number'
		 `new_date` date COMMENT 'Date entered into HGMD'
	*/
	if (starts_with($line, "INSERT INTO `grosins` VALUES"))
	{
		$line = substr($line, 31, -3);
		$entries = explode("'),('", $line);
		foreach($entries as $entry)
		{
			$entry = "'".$entry."'";
			//print $entry."\n";
			$parts = str_getcsv($entry, ',', '\'');
			//print_r($parts);
			if (count($parts)!=16) trigger_error("Expected 16 parts, but got ".count($parts)."!", E_USER_ERROR);
			
			list($disease, $gene, $type, $desc, $cdna, $tag, $author, $journal, $fullname, $vol, $page, $year, $pmid, $comments, $acc_num, $new_date) = $parts;
			
			$id = trim($acc_num);
			$class = trim($tag);
			$disease = trim($disease);
			$disease = strtr($disease, '"', '\'');
			$desc = trim($desc);
			$desc = strtr($desc, '"', '\'');
			$pmid = trim($pmid);
			if ($pmid=="NULL") $pmid = "";
			$type = trim($type);
			if ($type=="I") $type = "insertion";
			else if ($type=="D") $type = "duplication";
			else trigger_error("Invalid type '$type'!", E_USER_ERROR);
			
			//convert gene to coordinates
			list($chr, $start, $end, $gene_approved, $error) = gene2coords($gene);
			if (!is_null($error))
			{
				trigger_error($error." ($id)", E_USER_WARNING);
				continue;
			}
			print "{$chr}\t{$start}\t{$end}\t{$id} [GENE={$gene_approved} CLASS={$class} TYPE={$type} DISEASE=\"{$disease}\" DESC=\"{$desc}\" PMID={$pmid}]\n";
		}
	}
}
?>
