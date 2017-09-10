<?php

include("framework.php");

//##################################################################################
start_test("RTF::tablerows");

$result1 = rtf_table_rows(array(array("11","12"),array("21","22")),array(3000,6500),array("\qj","\qc"),true,false);
$result2 = rtf_table_rows(array(array("11","12"),array("21","22")),array(3000,6500),array("\qj","\qc"),false,true);

$expected1 = "\n{\trowd\trgaph70\fs20\clbrdrt\brdrw18\brdrs\clbrdrl\brdrw18\brdrs\clbrdrb\brdrw18\brdrs\clbrdrr\brdrw18\brdrs\cellx3000\clbrdrt\brdrw18\brdrs\clbrdrl\brdrw18\brdrs\clbrdrb\brdrw18\brdrs\clbrdrr\brdrw18\brdrs\cellx6500\pard\intbl\sa20\sb20\qj {11}\cell\pard\intbl\sa20\sb20\qc {12}\cell\row}\n{\trowd\trgaph70\fs20\clbrdrt\brdrw18\brdrs\clbrdrl\brdrw18\brdrs\clbrdrb\brdrw18\brdrs\clbrdrr\brdrw18\brdrs\cellx3000\clbrdrt\brdrw18\brdrs\clbrdrl\brdrw18\brdrs\clbrdrb\brdrw18\brdrs\clbrdrr\brdrw18\brdrs\cellx6500\pard\intbl\sa20\sb20\qj {21}\cell\pard\intbl\sa20\sb20\qc {22}\cell\row}";
check($result1, $expected1);

$expected2 = "\n{\trowd\trgaph70\fs20\cellx3000\cellx6500\pard\intbl\sa20\sb20\qj\b {11}\cell\pard\intbl\sa20\sb20\qc\b {12}\cell\row}\n{\trowd\trgaph70\fs20\cellx3000\cellx6500\pard\intbl\sa20\sb20\qj\b {21}\cell\pard\intbl\sa20\sb20\qc\b {22}\cell\row}";
check($result2, $expected2);
end_test();

//##################################################################################
start_test("RTF::tablerow");

$result1 = rtf_table_row(array("11","12"),array(3000,6500),array("\qj","\qc"),true,false);
$result2 = rtf_table_row(array("21","22"),array(3000,6500),array("\qj","\qc"),false,true);

$expected1 = "\n{\trowd\trgaph70\fs20\clbrdrt\brdrw18\brdrs\clbrdrl\brdrw18\brdrs\clbrdrb\brdrw18\brdrs\clbrdrr\brdrw18\brdrs\cellx3000\clbrdrt\brdrw18\brdrs\clbrdrl\brdrw18\brdrs\clbrdrb\brdrw18\brdrs\clbrdrr\brdrw18\brdrs\cellx6500\pard\intbl\sa20\sb20\qj {11}\cell\pard\intbl\sa20\sb20\qc {12}\cell\row}";
check($result1, $expected1);

$expected2 = "\n{\trowd\trgaph70\fs20\cellx3000\cellx6500\pard\intbl\sa20\sb20\qj\b {21}\cell\pard\intbl\sa20\sb20\qc\b {22}\cell\row}";
check($result2, $expected2);

end_test();

?>
