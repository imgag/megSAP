<?php

include("framework.php");

//##################################################################################
start_test("ToolBase::storeTDX");

//create dummy tool
$tool = new ToolBase("test_toolbase", "0.0a", "full description");
$tool->addInfile("in1", "desc1", false, false);
$tool->addInfile("in2", "desc2", true, false);
$tool->addOutfile("out1", "desc3", false, false);
$tool->addOutfile("out2", "desc4", true, false);
$tool->addInfileArray("in_array1", "desc5", false, false);
$tool->addInfileArray("in_array2", "desc6", true, false);
$tool->addFlag("flag1", "desc7");
$tool->addFlag("flag2", "desc8");
$tool->addString("string1", "desc9", false);
$tool->addString("string2", "desc10", true, "default_value");
$tool->addString("string3", "desc11", true);
$tool->addInt("int1", "desc12", false);
$tool->addInt("int2", "desc13", true, 4711);
$tool->addInt("int3", "desc14", true);
$tool->addFloat("float1", "desc15", false);
$tool->addFloat("float2", "desc16", true, 47.11);
$tool->addFloat("float3", "desc17", true);
$tool->addEnum("enum1", "desc18", false, array("a","b","c"));
$tool->addEnum("enum2", "desc19", true, array("d","e","f"), "d");
$tool->addEnum("enum3", "desc20", true, array("g","h","i"));

//store TDX file
$out_file = output_folder()."test_toolbase.tdx";
$tool->storeTDX($out_file);

//check well-formed
$xml_text = file_get_contents($out_file);
check(xml_is_wellformed($xml_text, $messages), true);

//check against schema
$xml_schema = file_get_contents(data_folder()."TDX_v1.xsd");
check(xml_matches_schema($xml_text, $xml_schema, $messages), true);

//extract version of ngs-bits
$ver = $tool->extractVersion(get_path("ngs-bits")."/SeqPurge");
check($ver!="n/a", true);

//exec
list($stdout, $stderr, $exit_code) = $tool->exec("echo", "bla && (echo bla2 >&2) && exit 17", false, false, false);
check(count($stdout)==1, true);
check(implode("", $stdout)=="bla", true);
check(count($stderr)==1, true);
check(implode("", $stderr)=="bla2", true);
check($exit_code==17, true);

end_test();

?>
