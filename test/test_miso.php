<?php

include("framework.php");

//##################################################################################
$miso = Miso::fromOBO();

start_test("MISO::TermByID");
$term_name = "transcript_variant";
$term_id = $miso->getTermIdByName($term_name);
check($term_id, "SO:0001576");
end_test();

start_test("MISO::ChildrenIDs");
$term_id = "SO:0001576";
$result = $miso->getTermChildrenIDs($term_id);
$expected = array("SO:0001621","SO:0001627","SO:0002011","SO:0001619","SO:0001791","SO:0001968","SO:0001577","SO:0001596","SO:0001568");
check(asort($result), asort($expected));
$term_id = "SO:0001627";
$result = $miso->getTermChildrenIDs($term_id,true);
$expected = array("SO:0001969","SO:0001629","SO:0002018","SO:0001970","SO:0001574","SO:0001575","SO:0001787");
check(asort($result), asort($expected));
end_test();

start_test("MISO::ChildrenNames");
$term_id = "SO:0001576";
$result = $miso->getTermChildrenNames($term_id);
$expected = array("SO:0001621","SO:0001627","SO:0002011","SO:0001619","SO:0001791","SO:0001968","SO:0001577","SO:0001596","SO:0001568");
check(asort($result), asort($expected));
$term_id = "SO:0001627";
$result = $miso->getTermChildrenNames($term_id,true);
$expected = array("SO:0001969","SO:0001629","SO:0002018","SO:0001970","SO:0001574","SO:0001575","SO:0001787");
check(asort($result), asort($expected));
end_test();

start_test("MISO::ValidTermByID");
$term_id = "SO:0001576";
$valid = Miso::isValidTermByID($term_id);
check($valid, true);
$term_id = "SO:0001579";
$valid = Miso::isValidTermByID($term_id);
check($valid, false);
end_test();

start_test("MISO::ValidTermByName");
$term_id = "SO:0001576";
$valid = Miso::isValidTermByName("transcript_variant");
check($valid, true);
$term_id = "SO:0001579";
$valid = Miso::isValidTermByName("test");
check($valid, false);
end_test();

?>
