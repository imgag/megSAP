<?php

include("framework.php");

//##################################################################################
start_test("ToolBase::storeTDX");

//create dummy tool
$tool = new ToolBase("test_toolbase", "0.0a", "full description");

//extract version
$ver = $tool->extractVersion("bgzip");
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
