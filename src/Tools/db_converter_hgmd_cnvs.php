<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$ngsbits = get_path("ngs-bits");

function gene2coords($gene)
{
	global $ngsbits;
	static $cache = [];
	
	if (!isset($cache[$gene]))
	{
		list($stdout, $stderr, $exit_code) = exec2("echo '{$gene}' | $ngsbits/GenesToApproved | cut -f1");
		$gene_approved = trim($stdout[0]);
		
		list($stdout, $stderr, $exit_code) = exec2("echo '{$gene_approved}' | $ngsbits/GenesToBed -source ensembl -mode gene | $ngsbits/BedMerge");
		$stdout = array_map('trim', $stdout);
		$stdout = array_diff($stdout, [""]);
		
		//error
		if (count($stdout)==0)
		{
			$cache[$gene] = [null, null, null, null, "No coordinates for gene '{$gene}'!"];
		}
		else if (count($stdout)>1)
		{
			trigger_error("Several genomic ranges for gene '{$gene}'. Selecting the largest range.", E_USER_WARNING);
			
			$max_size = -1;
			$max_size_i = -1;
			for($i=0; $i<count($stdout); ++$i)
			{
				list($chr, $start, $end) =  explode("\t", $stdout[0]);
				$size = $end-$start;
				if ($size>$max_size)
				{
					$max_size = $size;
					$max_size_i = $i;
				}
			}
			list($chr, $start, $end) =  explode("\t", $stdout[$max_size_i]);
			$cache[$gene] = [$chr, $start, $end, $gene_approved, null];
		}
		else
		{
			list($chr, $start, $end) =  explode("\t", $stdout[0]);
			$cache[$gene] = [$chr, $start, $end, $gene_approved, null];
		}
	}
	
	return $cache[$gene];
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
			if (count($parts)!=15) trigger_error("Expected 15 parts, but got ".count($parts)." in deletion entry: {$entry}", E_USER_ERROR);
			
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
			if (count($parts)!=16) trigger_error("Expected 16 parts, but got ".count($parts)." in insertion entry: {$entry}", E_USER_ERROR);
			
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
			else
			{
				trigger_error("Invalid type '$type' in insertion entry: {$entry}", E_USER_WARNING);
				continue;
			}
			
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
