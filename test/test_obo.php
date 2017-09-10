<?php

include("framework.php");

//##################################################################################
$file = repository_basedir()."/data/dbs/Ontologies/so-xp_3_0_0.obo";
$obo = Obo::fromOBO($file);

start_test("OBO::TermByID");
$term_name = "transcript_variant";
$term_id = $obo->getTermIdByName($term_name);
check($term_id, "SO:0001576");
end_test();

start_test("OBO::ChildrenIDs");
$term_id = "SO:0001576";
$result = $obo->getTermChildrenIDs($term_id);
$expected = array("SO:0001621","SO:0001627","SO:0002011","SO:0001619","SO:0001791","SO:0001968","SO:0001577","SO:0001596","SO:0001568");
check(asort($result), asort($expected));
$term_id = "SO:0001627";
$result = $obo->getTermChildrenIDs($term_id,true);
$expected = array("SO:0001969","SO:0001629","SO:0002018","SO:0001970","SO:0001574","SO:0001575","SO:0001787");
check(asort($result), asort($expected));
end_test();

start_test("OBO::ChildrenNames");
$term_id = "SO:0001576";
$result = $obo->getTermChildrenNames($term_id);
$expected = array("SO:0001621","SO:0001627","SO:0002011","SO:0001619","SO:0001791","SO:0001968","SO:0001577","SO:0001596","SO:0001568");
check(asort($result), asort($expected));
$term_id = "SO:0001627";
$result = $obo->getTermChildrenNames($term_id,true);
$expected = array("SO:0001969","SO:0001629","SO:0002018","SO:0001970","SO:0001574","SO:0001575","SO:0001787");
check(asort($result), asort($expected));
end_test();

start_test("OBO::ValidTermByID");
$term_id = "SO:0001576";
$valid = Obo::isValidTermByID($file, $term_id);
check($valid, true);
$term_id = "SO:0001579";
$valid = Obo::isValidTermByID($file, $term_id);
check($valid, false);
end_test();

start_test("OBO::ValidTermByName");
$term_id = "SO:0001576";
$valid = Obo::isValidTermByName($file, "transcript_variant");
check($valid, true);
$term_id = "SO:0001579";
$valid = Obo::isValidTermByName($file, "test");
check($valid, false);
end_test();

?>
