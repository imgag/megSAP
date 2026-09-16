<?php
declare(strict_types=1);
require_once(dirname($_SERVER['SCRIPT_FILENAME']) . "/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("convert_dragen_cram", "Converts a DRAGEN CRAM based on the Illumina GRCh38/pangenome into a CRAM compatible with the megSAP reference.");
$parser->addInfile("in", "Input CRAM file from DRAGEN.", false);
$parser->addInfile("in_ref", "Illumina linear reference genome FASTA which was used to create 'in'.", false);
$parser->addInfile("unsupported", "BED file with regions that are not present in megSAP reference genome. Reads in these regions will be unmapped in the output.", false);
$parser->addInt("threads", "Number of threads to use.", 1);
$parser->addOutfile("out", "Outut CRAM file for megSAP.", false);
$parser->addFlag("debug", "Enable debug output.");
extract($parser->parse($argv));

//extended desciption:
//- makes reads in unsupported regions unmapped and fixes mates if necessary
//- re-encode CRAM using megSAP genome
//  - renames chrM to chrMT
//  - fixes SAM tags where necessary

//init
$genome = genome_fasta("GRCh38");
$tmp_dir = $parser->tempFolder("dragen2megsap");

//extract contigs from a CRAM header
function parse_contigs_from_header($in_header): array
{
    $result = [];
    foreach ($in_header as $line)
    {
        if (!str_starts_with($line, "@SQ\t"))continue;

        $name = null;
        $length = null;

        foreach (explode("\t", $line) as $field)
        {
            if (str_starts_with($field, "SN:"))
            {
                $name = substr($field, 3);
            }
            elseif (str_starts_with($field, "LN:"))
            {
                $length = (int)substr($field, 3);
            }
        }

        if ($name!==null && $length!==null)
        {
            $result[$name] = $length;
        }
    }

    return $result;
}

//Convert DRAGEN contig names to megSAP names.
function to_megSAP_contig(string $chr): string
{
    if ($chr === 'chrM') return 'chrMT';
    return $chr;
}

//calculate reference span of a CIGAR string
function referenceSpan(string $cigar): int
{
    if ($cigar === '*') return 0;

    $length = 0;
    preg_match_all('/(\d+)([MIDNSHP=X])/', $cigar, $matches, PREG_SET_ORDER);
    foreach ($matches as $match)
    {
        if (in_array($match[2], ['M', 'D', 'N', '=', 'X'], true))
        {
            $length += (int)$match[1];
        }
    }

    return $length;
}

//checks if a SAM alignment overlaps a unsupported region
function alignmentIsUnsupported(string $chr, int $pos, string $cigar, array $unsupported): bool
{
    if ($chr === '*' || !isset($unsupported[$chr])) return false;

    $span = referenceSpan($cigar);
    if ($span === 0) return false;

    //SAM POS is 1-based; BED is 0-based
    $start = $pos - 1;
    $end = $start + $span;
	
	foreach($unsupported[$chr] as list($start2, $end2))
	{
		if (range_overlap($start, $end, $start2, $end2)) return true;
	}
	return false;
}

//process one QNAME group from the affected subset
function processAffectedGroup(array $group, array $unsupported, $keep_handle, $reset_handle, array &$stats): void
{
    foreach ($group as list($qname, $flags, $chr, $pos, $cigar, $line))
    {
		//discard all non-primary mappings (samtools fixmate does not always correct mate chromomes and then samtools sort throws an error...
        $mapped = ($flags & 0x4) === 0;
        $primary = ($flags & 0x900) === 0;
		if (!$mapped || !$primary)
		{
            ++$stats['discarded'];
            continue;
		}
		
        if (alignmentIsUnsupported($chr, $pos, $cigar, $unsupported))
        {
			fwrite($reset_handle, $line."\n");

			++$stats['reset'];

            continue;
        }
		
        fwrite($keep_handle, $line."\n");
		
        ++$stats['kept'];
    }
}

//check SA:Z: tag (locations of secondary/supplementary reads) for unsupported regions/chromosomes
function fix_SA_tag(string $tag, array $unsupported, array $target_contigs): ?string
{
    $entries = [];

    foreach (explode(';', substr($tag, 5)) as $entry)
    {
        if ($entry === '') continue;

        $fields = explode(',', $entry);
		
		//unknown/nonstandard SA entry: leave untouched
        if (count($fields) < 6)
        {
            $entries[] = $entry;
            continue;
        }

        $old_chr = $fields[0];

        if (alignmentIsUnsupported($old_chr, (int)$fields[1], $fields[3], $unsupported)) continue;

        $new_chr = to_megSAP_contig($old_chr);

        if (!isset($target_contigs[$new_chr])) continue;

        $fields[0] = $new_chr;

        $entries[] = implode(',', $fields);
    }

    if (count($entries) === 0) return null;

    return 'SA:Z:' . implode(';', $entries) . ';';
}

try
{	
	//make sure the Illumina genome FAI exists
    if (!file_exists(($in_ref.".fai")))
    {
		$parser->execApptainer("samtools", "samtools", "faidx $in_ref");
    }
	
	//get header of input CRAM
	$tmp_file = $parser->tempFile(".txt");
    $parser->execApptainer("samtools", "samtools", "view -H -T $in_ref $in > $tmp_file", [$in, $in_ref]);
	$in_header = [];
	foreach(file($tmp_file) as $line)
	{
		$line = nl_trim($line);
		if ($line=="") continue;
		
		$in_header[] = $line;
	}

    //get contigs of both worlds
    $input_contigs = parse_contigs_from_header($in_header);
    $target_contigs = genome_chr_sizes($genome);
	
	//-------------------------------------------------------------------------------
    //(1) create unsupported BED file
	//-------------------------------------------------------------------------------
	$unsupported_file = $unsupported;
    $unsupported = [];
    foreach(file($unsupported_file) as $line)
    {
        $line = nl_trim($line);
        if ($line=="" || starts_with($line, '#') || starts_with($line, 'track ') || starts_with($line, 'browser ')) continue;

        $parts = explode("\t", $line);
        if (count($parts) < 3) continue;

        list($chr, $start, $end) = $parts;

		//check contigs really exist
		if (!isset($input_contigs[$chr])) trigger_error("Contig '$chr' from unsupported BED not present in input CRAM!", E_USER_ERROR); 

        $unsupported[$chr][] = [$start, $end];
    }

    //check input contigs
    foreach ($input_contigs as $old_chr => $input_length)
    {
        $new_chr = to_megSAP_contig($old_chr);

        if (!isset($target_contigs[$new_chr]))
        {
			//trigger_error("Input contig '$old_chr' not supported by megSAP! Adding it to unsupported BED!", E_USER_WARNING);
			$unsupported[$old_chr][] = [0, $input_length];
			continue;
        }

        if ($target_contigs[$new_chr] !== $input_length) trigger_error("Contig '$old_chr' has differint length in input CRAM and megSAP reference!", E_USER_ERROR);
    }
	
	//merge unsupported BED
    $unsupported_regions = $tmp_dir.'/unsupported_regions.bed';
	$tmp_file = $parser->tempFile(".bed");
	$tmp = "";
	foreach($unsupported as $chr => $ranges)
	{
		foreach($ranges as list($start, $end))
		{
			$tmp .= "$chr\t$start\t$end\n";
		}
	}
	file_put_contents($tmp_file, $tmp);
	$parser->execApptainer("ngs-bits", "BedMerge", "-in $tmp_file -out $unsupported_regions");

	//-------------------------------------------------------------------------------
    //(2) extract read names
	//-------------------------------------------------------------------------------
	$time_start = microtime(true);
	print "extracting read names...\n";
	$unsupported_reads = $tmp_dir.'/unsupported_reads.sam';
	$parser->execApptainer("samtools", "samtools view", "-L $unsupported_regions -@ $threads  -F 0x4 -M -T $in_ref -o $unsupported_reads $in", [$unsupported_regions, $in_ref, $in]);
	$unsupported_read_ids = $tmp_dir.'/unsupported_read_ids.txt';
	exec2("cut -f1 {$unsupported_reads} | sort | uniq > {$unsupported_read_ids}");
	list($stdout) = exec2("wc -l < $unsupported_read_ids");
	$unsupported_read_count = trim(implode("", $stdout));
	if (!$debug) unlink($unsupported_reads);
	print "  extracted $unsupported_read_count read names in unsupported regions\n";
	print "  took ".time_readable(microtime(true)-$time_start)."\n";

	//-------------------------------------------------------------------------------
	//(3) extract affected reads by read names
	//-------------------------------------------------------------------------------
	$time_start = microtime(true);
	print "extracting affected reads from input CRAM...\n";
	$affected_bam = $tmp_dir.'/affected_reads.bam'; //reads overlapping unsupported BED
	$parser->execApptainer("samtools", "samtools view", "-@ $threads -T $in_ref -N $unsupported_read_ids -b -o $affected_bam $in", [$in, $in_ref]);
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
	//-------------------------------------------------------------------------------
	//(3) collate reads with the same name: they are next to each other in SAM then
	//-------------------------------------------------------------------------------
	$time_start = microtime(true);
	$affected_collate = $tmp_dir.'/affected_collate.sam';
	print "collating reads...\n";
	$parser->execApptainer("samtools", "samtools collate", "-@ $threads --output-fmt SAM -o $affected_collate $affected_bam");
	if (!$debug) unlink($affected_bam);
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
	//-------------------------------------------------------------------------------
	//(4) process collated reads - split into two file: keep vs. reset
	//-------------------------------------------------------------------------------
	print "processing collated reads...\n";
	$time_start = microtime(true);
	
	//open streams for tmp SAM files
	$keep_sam = $tmp_dir . '/affected.keep.sam';
	$reset_sam = $tmp_dir . '/affected.reset.sam';
	$keep_handle = fopen2($keep_sam, 'w');
	$reset_handle = fopen2($reset_sam, 'w');		

	//add header to tmp SAM files
	fwrite($keep_handle, implode("\n", $in_header)."\n");
	fwrite($reset_handle, implode("\n", $in_header)."\n");

	$stats = [
		'reset' => 0,
		'discarded' => 0,
		'kept' => 0
	];
	$current_qname = "";
	$group = [];
	$fp = fopen2($affected_collate, "r");
	while (!feof($fp))
	{
		$line = nl_trim(fgets($fp));
		if ($line=="" || $line[0]=="@") continue;
		
		$fields = explode("\t", nl_trim($line));
		if (count($fields) < 11) trigger_error("Invalid SAM record: $line", E_USER_ERROR);
		
		$qname = $fields[0];
		if ($qname!=$current_qname)
		{
			processAffectedGroup($group, $unsupported, $keep_handle, $reset_handle, $stats);
			$group = [];
		}

		$current_qname = $qname;
		$group[] = [$qname, (int)$fields[1], $fields[2], (int)$fields[3], $fields[5], $line];
	}
	processAffectedGroup($group, $unsupported, $keep_handle, $reset_handle, $stats);
	
	print "  count of primary alignments converted to unmapped: ".$stats['reset']."\n";
	print "  count of primary alignments kept: ".$stats['kept']."\n";
	print "  count of secondary/supplementary alignments discarded: ".$stats['discarded']."\n";
	
	//close file handles
	fclose($fp);
	fclose($keep_handle);
	fclose($reset_handle);
	if (!$debug) unlink($affected_collate);
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
	//-------------------------------------------------------------------------------
	//(5) Convert the two SAM streams back to BAM
	//-------------------------------------------------------------------------------
	print "converting SAM files to BAM...\n";
	$time_start = microtime(true);
	
	//convert keep SAM to BAM
	$keep_bam = $tmp_dir . '/affected.keep.bam';
	$parser->execApptainer("samtools", "samtools view", "-@ $threads -b -o $keep_bam $keep_sam");
	if (!$debug) unlink($keep_sam);
	
	//Correctly restores SEQ/QUAL orientation and removes alignment-specific state			
	$reset_bam = $tmp_dir . '/affected.reset.bam';
	$parser->execApptainer("samtools", "samtools reset", "-@ $threads --dupflag --output-fmt SAM -O BAM -o $reset_bam $reset_sam");
	if (!$debug) unlink($reset_sam);
	print "  took ".time_readable(microtime(true)-$time_start)."\n";

	//-------------------------------------------------------------------------------
	//(6) Recombine affected subset, then repair mate fields
	//-------------------------------------------------------------------------------
	print "re-combining affected read BAMs...\n";
	$time_start = microtime(true);
	
	//merge
	$affected_remerged = $tmp_dir . '/affected.remerged.bam';
	$parser->execApptainer("samtools", "samtools cat", "$keep_bam $reset_bam -o $affected_remerged");
	if (!$debug) unlink($keep_bam);
	if (!$debug) unlink($reset_bam);
	
	//collate again (needed as we separately handled reset/keeep reads before)
	$affected_recollated = $tmp_dir . '/affected.recollated.bam';
	$parser->execApptainer("samtools", "samtools collate", "-@ $threads -o $affected_recollated $affected_remerged");
	if (!$debug) unlink($affected_remerged);
	
	//fix mates infos
	$affected_fixed = $tmp_dir . '/affected.fixed.bam';
	$parser->execApptainer("samtools", "samtools fixmate", "-@ $threads $affected_recollated $affected_fixed");
	if (!$debug) unlink($affected_recollated);

	//sort (note: we cannot used ToolBase::sortBam because we need to use a different reference genome)
	$affected_resorted = $tmp_dir . '/affected.resorted.bam';
	$parser->execApptainer("samtools", "samtools sort", "-T $tmp_dir/sort_fixed -@ {$threads} -m 1G -T $in_ref -o $affected_resorted $affected_fixed", [$in_ref]);
	if (!$debug) unlink($affected_fixed);
	print "  took ".time_readable(microtime(true)-$time_start)."\n";

	//-------------------------------------------------------------------------------
	//(7) re-combine unaffected and affected reads
	//-------------------------------------------------------------------------------
	print "re-combining affected und unaffected read BAMs...\n";
	$time_start = microtime(true);
	
	//extract reads in supported regions
	$unaffected_bam = $tmp_dir . '/unaffected.bam';
	$parser->execApptainer("samtools", "samtools view", "-@ {$threads} -T $in_ref -N ^$unsupported_read_ids -o $unaffected_bam $in", [$in_ref, $in]);
	
	//merge unaffected and affected reads
	$combined_bam = $tmp_dir . '/combined.bam';
	$parser->execApptainer("samtools", "samtools merge", "-c -@ $threads -f -h $unaffected_bam --reference $in_ref -o $combined_bam $affected_resorted $unaffected_bam", [$in_ref, $affected_resorted, $unaffected_bam]);
	if (!$debug) unlink($affected_resorted);
	if (!$debug) unlink($unaffected_bam);
	print "  took ".time_readable(microtime(true)-$time_start)."\n";

	//-------------------------------------------------------------------------------
	//(8) Build target header
	//-------------------------------------------------------------------------------
	print "building target header...\n";
	$time_start = microtime(true);
	
	//add format header
    $target_header = [];
    foreach ($in_header as $line)
    {
        if (str_starts_with($line, '@HD'))
        {
            $target_header[] = $line;
            break;
        }
    }
    if (count($target_header) === 0) trigger_error("No @HD line input CRAM file!", E_USER_ERROR);

    //use megSAP @SQ dictionary
	$target_dict = $tmp_dir . '/target_dict.txt';
    $parser->execApptainer("samtools", "samtools dict > $target_dict", $genome, [$genome]);
    foreach (file($target_dict) as $line)
    {
		$line = nl_trim($line);
        if (str_starts_with($line, '@SQ'))
        {
            $target_header[] = $line;
        }
    }

    //preserve other header data from DRAGEN input CRAM
    foreach ($in_header as $line)
    {
        if ($line === '' || str_starts_with($line, '@HD') || str_starts_with($line, '@SQ')) continue;

        $target_header[] = $line;
    }
	
    $target_header[] = "@CO\tConverted from DRAGEN reference to megSAP reference; unsupported primary alignments converted to unmapped!";
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
	//-------------------------------------------------------------------------------
	//(9) Rewrite alignment records to SAM using the new dictionary
	//-------------------------------------------------------------------------------
	print "rewriting BAM based on target reference (rename chrM and fix SAM tags) and sort it...\n";
	$time_start = microtime(true);
	
	//open output process
	$half_threads = round($threads/2);
	$command_sort = $parser->execApptainer("samtools", "samtools sort", "-@ $half_threads -m 1G -T $tmp_dir/sort_combined --reference $genome -u -", [$genome], [], true);
	$command_cram = $parser->execApptainer("samtools", "samtools view", "-@ $half_threads -C --output-fmt-option version=3.0 -T $genome -o $out -", [$genome], [dirname($out)], true);
	
	$p_out = popen("bash -c \"set -o pipefail && $command_sort | $command_cram\"", 'w');	
    if ($p_out===false) trigger_error("Could not start 'samtools sort | samtools view' output process!", E_USER_ERROR);
	
	//write header
	foreach ($target_header as $line)
	{
		fwrite($p_out, $line . "\n");
	}
	
	//read BAM with input process
    $p_in = popen("samtools view $combined_bam", 'r');
    if ($p_in===false) trigger_error("Could not start 'samtools view' reading process on BAM!", E_USER_ERROR);
    while (($line = fgets($p_in)) !== false)
    {
        $line = nl_trim($line);
        $fields = explode("\t", $line);
        if (count($fields) < 11) trigger_error("Invalid SAM record: $line", E_USER_ERROR);

		//fix read contig
        $chr = $fields[2];
		if ($chr!=='*')
        {
            $chr = to_megSAP_contig($chr);
            if (!isset($target_contigs[$chr])) trigger_error("Invalid read contig $chr", E_USER_ERROR);
			$fields[2] = $chr;
        }

		//fix mate contig
		$chr_mate = $fields[6];
        if ($chr_mate!== '*' && $chr_mate!=='=')
        {
            $chr_mate = to_megSAP_contig($chr_mate);
            if (!isset($target_contigs[$chr_mate])) trigger_error("Invalid mate contig $chr_mate", E_USER_ERROR);
			$fields[6] = $chr_mate;
        }

        //fix tags
        $tags = [];
        for ($i = 11; $i < count($fields); ++$i)
        {
            $tag = $fields[$i];

            //rename chrM -> chrMT and remove references to deleted alignments.
            if (str_starts_with($tag, 'SA:Z:'))
            {
                $tag = fix_SA_tag($tag, $unsupported, $target_contigs);
                if ($tag===null) continue;
            }

            //these tags can contain coordinates/reference names from the old reference. They are not needed by megSAP and are safer to remove.
            if (str_starts_with($tag, 'XA:Z:') || str_starts_with($tag, 'OA:Z:') || str_starts_with($tag, 'CC:Z:') || str_starts_with($tag, 'CP:i:')) continue;
          
            $tags[] = $tag;
        }

        $fields = array_merge(array_slice($fields, 0, 11), $tags);
        fwrite($p_out, implode("\t", $fields)."\n");
    }
    $exit_code = pclose($p_in);
	if ($exit_code !== 0) trigger_error("'samtools view' reading process failed with exit code '$exit_code'", E_USER_ERROR);
    $exit_code = pclose($p_out);
	if ($exit_code !== 0) trigger_error("'samtools sort | samtools view' writing process failed with exit code '$exit_code'", E_USER_ERROR);
	if (!$debug) unlink($combined_bam);
	
	//index
	$parser->indexBam($out, $threads);
	
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
}
finally
{
    exec2("rm -rf $tmp_dir");
}

?>