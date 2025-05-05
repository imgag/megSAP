<?php
/** 
	@page mvh_grz_export
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("mvh_grz_export", "GRZ export for Modellvorhaben.");
$parser->addInt("case_id", "Case ID in Modellvorhaben DB.", false);
extract($parser->parse($argv));

//init
$db_mvh = DB::getInstance("MVH");
$db_ngsd = DB::getInstance("NGSD");
$mvh_folder = get_path("mvh_folder");

//check that case ID is valid
$seq_type = $db_mvh->getValue("SELECT seq_type FROM case_data WHERE id='{$case_id}'", "");
if ($seq_type=="") trigger_error("No case with id '{$case}' in MVH database!", E_USER_ERROR);
if ($seq_type!="wes") trigger_error("Type '{$seq_type}' not supported yet!", E_USER_ERROR); //TODO implement WGS, lrGS, TN WES, TN WGS
$disease_type = $db_mvh->getValue("SELECT disease_type FROM case_data WHERE id='{$case_id}'");
print "case: {$case_id} (seq_type: {$seq_type}, disease_type: {$disease_type})\n";

//determine KDK
$kdk = "";
if ($disease_type=="rare") $kdk = "KDKTUE002"; //NSE - Tübingen
if ($disease_type=="oncological") $kdk = "KDKTUE005"; //DNPM - Tübingen
if ($disease_type=="hboc") $kdk = "KDKL00003"; //DK-FBREK - Leipzig
if ($disease_type=="hereditary_cancer") $kdk = "KDKL00004"; //DK-FDK - Leipzig
if ($kdk=="") trigger_error("Could not determine KDK for sample '{$ps1}' with disease_type '{$disease_type}'!", E_USER_ERROR);
	
//check germline processed sample is ok
$ps1 = $db_mvh->getValue("SELECT ps1 FROM case_data WHERE id='{$case_id}'");
$info1 = get_processed_sample_info($db_ngsd, $ps1);
$sys1 = $info1['sys_name_short'];
$patient_id = $info1['patient_identifier'];
if ($patient_id=="") trigger_error("No patient identifier set for sample '{$ps1}''!", E_USER_ERROR);
print "sample 1: {$ps1} (system: {$sys1})\n";
$readl1 = get_processed_sample_qc($db_ngsd, $ps1, "QC:2000006");
while (contains($readl1, "-")) $readl1 = explode("-", $readl1, 2)[1];

//determine ROI
$roi = "";
if ($sys1=="twistCustomExomeV2" || $sys1=="twistCustomExomeV2Covaris") $roi = "{$mvh_folder}/rois/twist_exome_core_plus_refseq.bed";
if ($sys1=="hpHBOCv5") $roi = "/mnt/storage2/megSAP/data/enrichment/hpHBOCv5_2014_10_27.bed"; //TODO: only for testing > remove!
if (!contains($seq_type, "wgs") && !contains($seq_type,"lrgs") && $roi=="") trigger_error("Could not determine target region for sample '{$ps1}'' with processing system '{$sys1}'!", E_USER_ERROR);

//determine study sub-type
$study_subtype = "";
if ($seq_type=="wes" || $seq_type=="wgs" || $seq_type=="lrgs") $study_subtype = "germline-only";
if ($seq_type=="wes_tn" || $seq_type=="wgs_tn") $study_subtype = "tumor+germline";
if ($study_subtype=="") trigger_error("Could not determine strudy subtype for case {$case_id}!", E_USER_ERROR);

//check tumor sample
//TODO exists
//TODO system same as germline sample

//determine tanG==VNg
$sub_ids = $db_mvh->getValues("SELECT id FROM `submission_grz` WHERE status='pending'");
if (count($sub_ids)!=1)  trigger_error(count($sub_ids)." pending GRZ submissions for case {$case_id}. Must be one!", E_USER_ERROR);
$sub_id = $sub_ids[0];
$tang = $db_mvh->getValue("SELECT tang FROM submission_grz WHERE id='{$sub_id}'");

//create export folder
$folder = "{$mvh_folder}/export_grz/{$case_id}/";
exec2("mkdir -p {$folder}");
print "export folder: {$folder}\n";
exec2("mkdir -p {$folder}/files/");
exec2("mkdir -p {$folder}/metadata/");

//create raw data (FASTQs + target region)
$bam1 = $info1['ps_bam'];
$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$bam1} -out1 {$folder}/files/{$tang}.read1.fastq.gz -out2 {$folder}/files/{$tang}.read2.fastq.gz -extend {$readl1}", [$bam1], ["{$folder}/files/"]);
if ($roi!="") exec2("cp {$roi} {$folder}/files/target_regions.bed");

//create meta data
$json = [];
$json[] = "{\n";
$json[] = "\"\$schema\": \"https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/v1.1.7/GRZ/grz-schema.json\",\n";
$json[] = "  \"submission\": {\n";
$json[] = "    \"submissionDate\": \"".$db_mvh->getValue("SELECT date FROM submission_grz WHERE id='{$sub_id}'")."\",\n";
$json[] = "    \"submissionType\": \"".$db_mvh->getValue("SELECT type FROM submission_grz WHERE id='{$sub_id}'")."\",\n";
$json[] = "    \"submitterId\": \"260840108\",\n";
$json[] = "    \"tanG\": \"{$tang}\",\n";
$json[] = "    \"localCaseId\": \"{$patient_id}\",\n";
$json[] = "    \"genomicDataCenterId\": \"GRZTUE002\",\n";
$json[] = "    \"clinicalDataNodeId\": \"{$kdk}\",\n";
$json[] = "    \"labName\": \"Institute of Medical Genetics and Applied Genomics, Tuebingen, Germany\",\n";
$json[] = "    \"genomicStudyType\": \"single\",\n";
$json[] = "    \"genomicStudySubtype\": \"{$study_subtype}\",\n";
$json[] = "    \"coverageType\": \"UNK\",\n"; //TODO get info and implement: "GKV" --> Gesetzliche Krankenversicherung "PKV" --> Private Krankenversicherung "BG", --> Berufsgenossenschaft "SEL", --> Selbstzahler "SOZ", --> Sozialamt "GPV", --> Gesetzliche Pflegeversicherung "PPV", --> Private Pflegeversicherung "BEI", --> Beihilfe "SKT", --> Sonstige Kostenträger "UNK" --> Unklar/Unbekannt
$json[] = "    \"diseaseType\": \"{$disease_type}\"\n";
$json[] = "  },\n";
$json[] = "  \"donors\":\n";
$json[] = "  [\n";
$json[] = "  ]\n";
$json[] = "{\n";
print_r($json);


//Validate the submission
#grz-cli validate --submission-dir EXAMPLE_SUBMISSION

//Encrypt the submission
#grz-cli encrypt --submission-dir EXAMPLE_SUBMISSION

//Upload the submission
#grz-cli upload --submission-dir EXAMPLE_SUBMISSION


/*
TODO:
- test if export works with other user, e.g. bioinf
- store schema of MVH DB


#Installation notes GRZ-CLI (see https://github.com/BfArM-MVH/grz-cli):

	- Install miniforge at /mnt/storage2/megSAP/tools/miniforge3/
		> curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
		> bash Miniforge3-$(uname)-$(uname -m).sh
	- Install GRZ-CLI
		> /mnt/storage2/megSAP/tools/miniforge3/bin/conda create -n grz-tools -c conda-forge -c bioconda "grz-cli"
		> /mnt/storage2/megSAP/tools/miniforge3/bin/conda activate grz-tools
	- Updates with:
		> /mnt/storage2/megSAP/tools/miniforge3/bin/conda update -n grz-tools "grz-cli"
	- Activate:
		> /mnt/storage2/megSAP/tools/miniforge3/bin/conda activate grz-tools
*/
?>
