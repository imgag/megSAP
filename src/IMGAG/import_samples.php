<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("import_samples", "Import samples received through https://megsap.de/ngs_data_transfer/");
$parser->addInfileArray("xml", "XML file(s) for same project/sender.", false);
$parser->addFlag("check_only", "Perform XML checks only. Do not import into NGSD.");

extract($parser->parse($argv));

//helper function for XML booleans
function xmlBool(string $value): bool
{
    return match ($value)
    {
        'true', '1' => true,
        'false', '0' => false,
        default => throw new RuntimeException("Invalid XML boolean value: '$value'"),
    };
}

//helper function to parse XML into associative array
function parseNgsMetaData(string $xml_file): array
{
    libxml_use_internal_errors(true);

    $dom = new DOMDocument();

    if (!$dom->load($xml_file, LIBXML_NONET))
    {
        $errors = libxml_get_errors();
        libxml_clear_errors();

        throw new RuntimeException(
            "Could not parse XML:\n" .
            implode(
                "\n",
                array_map(
                    fn($e) => trim($e->message),
                    $errors
                )
            )
        );
    }

    $xml = simplexml_import_dom($dom);

    if ($xml === false)
    {
        throw new RuntimeException(
            "Could not convert DOM to SimpleXML"
        );
    }


    $result = [
        'version' => (string)$xml['version'],

        'sender' => [
            'name' => (string)$xml->Sender['name'],
            'affiliation' => (string)$xml->Sender['affiliation'],
            'email' => (string)$xml->Sender['email'],
        ],

        'project' => [
            'name' => (string)$xml->Project['name'],
        ],

        'samples' => [],
    ];


    foreach ($xml->Sample as $sample)
    {
        $patient_info = $sample->PatientInfo;

        $patient = [
            'gender' => (string)$patient_info['gender'],
            'disease_group' => (string)$patient_info['disease_group'],
            'disease_status' => (string)$patient_info['disease_status'],
            'yob' => isset($patient_info['yob'])
                ? (int)$patient_info['yob']
                : null,
            'hpo' => [],
        ];


        foreach ($patient_info->HPO as $hpo)
        {
            $patient['hpo'][] = (string)$hpo['id'];
        }


        $files = [];

        foreach ($sample->File as $file)
        {
            $files[] = [
                'name' => (string)$file['name'],
                'type' => (string)$file['type'],
            ];
        }


        $relations = [];

        foreach ($sample->Relation as $relation)
        {
            $relations[] = [
                'relation' => (string)$relation['relation'],
                'sample_id' => (string)$relation['sample_id'],
            ];
        }


        $entry = [
            'id' => (string)$sample['id'],
            'type' => (string)$sample['type'],
            'tissue' => (string)$sample['tissue'],

            'is_tumor' => xmlBool(
                (string)$sample['is_tumor']
            ),

            'is_ffpe' => xmlBool(
                (string)$sample['is_ffpe']
            ),

            'sequencing_platform' =>
                (string)$sample['sequencing_platform'],

            'sequencing_type' =>
                (string)$sample['sequencing_type'],

            'sequencer' =>
                (string)$sample['sequencer'],

            'kit' =>
                (string)$sample['kit'],

            'adapter1_p5' =>
                isset($sample['adapter1_p5'])
                    ? (string)$sample['adapter1_p5']
                    : null,

            'adapter2_p7' =>
                isset($sample['adapter2_p7'])
                    ? (string)$sample['adapter2_p7']
                    : null,

            'patient' => $patient,
            'files' => $files,
            'relations' => $relations,
        ];


        $result['samples'][] = $entry;
    }


    return $result;
}

//get checksum from a file created by sha256sum
function get_checksum(string $filename): string
{
	$checksum = "";
	foreach(file($filename) as $line)
	{
		$line = str_simplified($line);
		if ($line=="") continue;
		if ($checksum!="") trigger_error("Checksum file '$filename' contains more than one checksum!", E_USER_ERROR);
		$checksum = explode(" ", $line)[0];
	}
	
	if ($checksum=="") trigger_error("Checksum file '$filename' contains no checksum!", E_USER_ERROR);
	return $checksum;
}

//get next external sample name
function next_external_sample_name(): string
{
	global $db;
	
	$prefix = "EX".date('y');
	$snext = 1;
	foreach($db->getValues("SELECT name FROM sample WHERE name LIKE '{$prefix}%' ORDER BY name ASC") as $sname)
	{
		$sname = substr($sname, strlen($prefix));
		$sname = ltrim($sname, '0');
		$snext = max($sname+1, $snext);
	}
	$snext = $prefix.str_pad($snext, 4, '0',STR_PAD_LEFT);
	
	return $snext;
}

//init
$schema = $parser->tempFile(".xsd");
exec2("wget -O - https://raw.githubusercontent.com/imgag/NGS_data_transfer/refs/heads/main/NgsMetaData_schema.xsd > $schema");
$db = DB::getInstance("NGSD");

//parse meta data from xml
print "Parsing XML file(s)...\n";
$project = null;
$sender = null;
$samples = [];
foreach($xml as $filename)
{
	//check XML is valid
	list($stdout, $stderr, $exit_code) = exec2("xmllint --noout --schema $schema $filename 2>&1", false);
	if ($exit_code!=0)
	{
		trigger_error("XML file $filename is not valid:\n".implode("\n", $stdout), E_USER_ERROR);
	}
	
	$data = parseNgsMetaData($filename);
	
	//project
	$tmp = $data['project']['name'];
	if ($project==null)
	{
		$project = $tmp;
	}
	else if ($tmp!=$project) trigger_error("Mismatching projects!", E_USER_ERROR);
	
	//sender
	$tmp = [$data['sender']['name'], $data['sender']['affiliation'], $data['sender']['email']];
	if ($sender==null)
	{
		$sender = $tmp;
	}
	else if ($tmp!=$sender) trigger_error("Mismatching senders!", E_USER_ERROR);
	
	//samples
	foreach($data['samples'] as $tmp)
	{
		$samples[] = $tmp;
	}
}

//check that sample relations are valid
print "Checking sample relations...\n";
$sample_ids = [];
foreach($samples as $sample)
{
	$sample_ids[$sample['id']] = true;
}
foreach($samples as $sample)
{
	foreach($sample['relations'] as $relation)
	{
		$id = $relation['sample_id'];
		if (!isset($sample_ids[$id])) trigger_error("Referenced sample '$id' not in sample list!", E_USER_ERROR);
	}
}

//check checksums
print "validating checksums...\n";
$folder = dirname(realpath($xml[0]))."/";
foreach($samples as $sample)
{
	foreach($sample['files'] as $file)
	{
		if ($file['type']!="BAM") trigger_error("Unsupported file type '".$file['type']."'!", E_USER_ERROR);
		
		//check file and external checksum are there
		$filename = $folder.$file['name'];
		if (!file_exists($filename)) trigger_error("File '$filename' not found!", E_USER_ERROR);
		$checksum1 = $filename.".sha256";
		if (!file_exists($checksum1)) trigger_error("Checksum file '$checksum1' not found!", E_USER_ERROR);
		
		//re-create checksum 
		$checksum2 = $filename.".sha256.imgag";
		if (!file_exists($checksum2))
		{
			print "  creating checksum of $filename\n";
			exec2("sha256sum $filename > $checksum2");
		}
		
		//validate checksums match
		if (get_checksum($checksum1)!=get_checksum($checksum2)) trigger_error("Checksum mismatch for '$filename'!", E_USER_ERROR);
	}
}

//get user id
$user_name = trim(exec('whoami'));
$user_id = $db->getValue("SELECT id FROM user WHERE user_id='$user_name'", -1);
if ($user_id==-1) trigger_error("User '$user_name' not in NGSD!", E_USER_ERROR);
print "user_id: $user_id\n";

//get project ID from NGSD
$project_id = $db->getValue("SELECT id FROM project WHERE name='$project'", -1);
if ($project_id==-1) trigger_error("Project '$project' not in NGSD!", E_USER_ERROR);
print "project_id: $project_id\n";

//get sender ID from NGSD (insert sender if necessary)
$sender_id = $db->getValue("SELECT id FROM sender WHERE name='".$sender[0]."'", -1);
if ($sender_id==-1)
{
	$db->executeStmt("INSERT INTO `sender` (`name`, `affiliation`, `email`) VALUES ('".implode("', '", $sender)."')");
	$sender_id = $db->getValue("SELECT id FROM sender WHERE name='".$sender[0]."'", -1);
}
print "sender_id: $sender_id\n";

//get species id
$species_id = $db->getValue("SELECT id FROM `species` WHERE `name` LIKE 'Human'");
print "species_id: $species_id\n";

//get run id
$run_id = $db->getValue("SELECT id FROM `sequencing_run` WHERE `name` LIKE '#00000'");
print "run_id: $run_id\n";

//check only
if ($check_only)
{
	print "Skipping NGSD import because of flag 'check_only'.\n";
	exit(0);
}

//import sample meta data
$date = date('Y-m-d');
$sample_ids = [];
$sampleid2ngsd = [];
foreach($samples as $sample)
{
	print "Importing sample ".$sample['id']."...\n";

	//get/create processing system
	$type = $sample['sequencing_type'];
	$type = strtr($type, ["srWGS"=>"WGS", "lrWGS"=>"lrGS"]);
	$sys_id = $db->getValue("SELECT id FROM `processing_system` WHERE `name_manufacturer` LIKE '".$sample['kit']."' AND `type` LIKE '".$type."' AND `platform` LIKE '".$sample['sequencing_platform']."'", -1);
	if ($sys_id==-1)
	{
		print "  Adding new processing system for kit '".$sample['kit']."'\n";
		exit(1); //TODO add including adapters
	}
	
	$name = next_external_sample_name();
	
	//insert sample
	$fields = [];
	$values = [];
	$fields[] = "name";
	$values[] = $name;
	$fields[] = "name_external";
	$values[] = $sample['id'];
	$fields[] = "received";
	$values[] = $date;
	$fields[] = "receiver_id";
	$values[] = $user_id;
	$fields[] = "sample_type";
	$values[] = $sample['type'];
	$fields[] = "tissue";
	$values[] = $sample['tissue'];
	$fields[] = "species_id";
	$values[] = $species_id;
	$fields[] = "gender";
	$values[] = $sample['patient']['gender'];
	$fields[] = "tumor";
	$values[] = $sample['is_tumor'] ? "1" : "0";
	$fields[] = "ffpe";
	$values[] = $sample['is_ffpe'] ? "1" : "0";
	$fields[] = "sender_id";
	$values[] = $sender_id;
	$fields[] = "disease_group";
	$values[] = $sample['patient']['disease_group'];
	$fields[] = "disease_status";
	$values[] = $sample['patient']['disease_status'];
	if ($sample['patient']['yob']!="")
	{
		$fields[] = "year_of_birth";
		$values[] = $sample['patient']['yob'];
	}
	$fields[] = "comment";
	$values[] = "Imported using import_samples.php by $user_name on $date";
	$db->executeQuery("INSERT INTO `sample` (`".implode("`, `", $fields)."`) VALUES ('".implode("', '", $values)."')");
	$sample_id = $db->getValue("SELECT id FROM sample WHERE name LIKE '$name'");
	print "  name: $name\n";
	print "  sample_id: $sample_id\n";
	$sample_ids[] = $sample_id;
	$sampleid2ngsd[$sample['id']] = [$sample_id, $sys_id, $name];
	
	//add HPO terms
	foreach($sample['patient']['hpo'] as $hpo_id)
	{
		$id = $db->getValue("SELECT id from hpo_term WHERE hpo_id='$hpo_id'", -1);
		if ($id==-1)
		{
			print "  Could not import unknown HPO term with ID '$hpo_id'\n";
			continue;
		}
		$db->executeQuery("INSERT INTO `sample_disease_info`(`sample_id`, `disease_info`, `type`, `user_id`, `date`) VALUES ($sample_id,'$hpo_id','HPO term id',$user_id,'$date')");
	}
}

//for debugging
/*
$tmp = implode(", ", $sample_ids);
print "TO DELETE:\n";
print "  DELETE FROM sample_disease_info WHERE sample_id IN ($tmp);\n";
print "  DELETE FROM sample_relations WHERE sample1_id IN ($tmp);\n";
print "  DELETE FROM sample_relations WHERE sample2_id IN ($tmp);\n";
print "  DELETE FROM processed_sample WHERE sample_id IN ($tmp);\n";
print "  DELETE FROM sample WHERE id IN ($tmp);\n";
*/

//add sample relations
print "Adding sample relations...\n";
foreach($samples as $sample)
{
	$sample_id = $sample['id'];
	$s1 = $sampleid2ngsd[$sample_id][0];
	foreach($sample['relations'] as $relation)
	{
		$related_id = $relation['sample_id'];
		$s2 = $sampleid2ngsd[$related_id][0];
		
		$relation = $relation['relation'];
		if ($relation=="sibling of")
		{
			$relation2 = "siblings";
		}
		else if ($relation=="child of")
		{
			$relation2 = "parent-child";
		}
		else trigger_error("Unhandled relation '$relation'!", E_USER_ERROR);
		
		$db->executeQuery("INSERT INTO `sample_relations`(`sample1_id`, `relation`, `sample2_id`, `user_id`) VALUES ('$s2','$relation2','$s1','$user_id')");
	}
}

//add processed samples
print "Adding processed samples...\n";
foreach($samples as $sample)
{
	$sample_id = $sample['id'];
	
	$fields = [];
	$values = [];
	$fields[] = "sample_id";
	$values[] = $sampleid2ngsd[$sample_id][0];
	$fields[] = "process_id";
	$values[] = "1";
	$fields[] = "sequencing_run_id";
	$values[] = $run_id;
	$fields[] = "lane";
	$values[] = "1";
	$fields[] = "processing_system_id";
	$values[] = $sampleid2ngsd[$sample_id][1];
	$fields[] = "project_id";
	$values[] = $project_id;	
	$db->executeQuery("INSERT INTO `processed_sample` (`".implode("`, `", $fields)."`) VALUES ('".implode("', '", $values)."')");
}

//move BAM files to sample folder
foreach($samples as $sample)
{
	$sample_id = $sample['id'];
	print "Copying BAM file of sample $sample_id ...\n";
		
	//create folder
	$ps = $sampleid2ngsd[$sample_id][2]."_01";
	list($stdout) = exec2("SamplePath -ps {$ps} -type SAMPLE_FOLDER");
	$ps_folder = trim(implode("", $stdout));
	exec2("mkdir -p $ps_folder", false);
	exec2("chmod -R 777 $ps_folder", false);
	
	foreach($sample['files'] as $file)
	{
		if ($file['type']!="BAM") trigger_error("Unsupported file type '".$file['type']."'!", E_USER_ERROR);

		$filename = $folder.$file['name'];
		$bam = $ps_folder."{$ps}.bam";
		exec2("cp $filename $bam");
		exec2("chmod 777 $bam", false);
	}
}

?>