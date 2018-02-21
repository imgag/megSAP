<?php
/** 
	@page vcf_merge
	
	@todo merge without vcftools 'vcf-merge'
	@todo make format/sample column content configurable
	@todo make quality/info column content configurable
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("vcf_merge", "Merges VCF files to one multi-sample VCF file.");
$parser->addInfileArray("in", "Input VCF.GZ files.", false);
$parser->addOutfile("out", "Output VCF file.", false);
extract($parser->parse($argv));

//create tmp folder
$tmp_folder = temp_folder("vcf_merge");
print "Temporary folder for simlyfied VCFs: $tmp_folder\n";

//create simplyfied files
$simple =  array();
foreach($in as $filename)
{
	print "Simplyfying $filename...\n";
	if (!ends_with($filename, ".vcf.gz"))
	{
		trigger_error("Input file $filename does have '.vcf.gz' extension.", E_USER_ERROR);
	}
	$output = array();
	list($file) = exec2("zcat $filename");
	foreach($file as $line)
	{
		$line = trim($line);
		if($line=="") continue;
		
		//header
		if($line[0]=="#")
		{
			$output[] = $line."\n";
			continue;
		}
		
		//content
		list($chr, $pos, $id, $ref, $obs, $qual, $filter, $info, $format, $sample) = explode("\t", $line);
		
		//fix info (remove DP)
		$info = explode(";", $info);
		$info2 = array();
		for($i=0; $i<count($info); ++$i)
		{
			if (starts_with($info[$i], "DP=")) continue;
			$info2[] = $info[$i];
		}
		$info = implode(";", $info2);
		
		//extract format info (GT)
		$format = explode(":", $format);
		$sample = explode(":", $sample);
		$format2 = array();
		$sample2 = array();
		for($i=0; $i<count($format); ++$i)
		{
			if ($format[$i]=="GT")
			{
				$sample[$i] = strtr($sample[$i], array("."=>"0", "|"=>"/"));
				$format2[] = $format[$i];
				$sample2[] = $sample[$i];
			}
		}
		$format = implode(":", $format2);
		$sample = implode(":", $sample2);
		$output[] = "$chr\t$pos\t$id\t$ref\t$obs\t$qual\t$filter\t$info\t$format\t$sample\n";
	}

	//store output
	$out_vcf = $tmp_folder."/".basename($filename, ".gz");
	file_put_contents($out_vcf, $output);
	
	//bgzip + index
	exec2("bgzip $out_vcf");
	$out_vcfgz = $out_vcf.".gz";
	exec2("tabix -p vcf $out_vcfgz");
	$simple[] = $out_vcfgz;
}


//merge simplyfied files
print "Merging...\n";
putenv("PERL5LIB=/mnt/share/opt/vcftools-78add55-bin/share/perl/5.14.2/:".getenv("PERL5LIB"));
exec2("/mnt/share/opt/vcftools-78add55-bin/bin/vcf-merge --ref-for-missing 0/0 ".implode(" ", $simple)." > $out");

//cleanup
exec2("rm -rf $tmp_folder");

?>