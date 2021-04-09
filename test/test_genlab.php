<?php

include("framework.php");


start_test("GenLabDB");


//skip tests if not enabled
if (GenLabDB::isEnabled())
{
	//create
	$db = GenLabDB::getInstance();
	
	//quote
	check($db->quote("'"), "''''");
	
	//executeQuery	
	$res = $db->executeQuery("SELECT TOP 10 * FROM v_ngs_duo ORDER BY Labornummer_Index");
	check(count($res), 10);
	check(count($res[0]), 3);
	
	//getValue
	$res = $db->getValue("SELECT TOP 1 Labornummer_Index FROM v_ngs_duo ORDER BY Labornummer_Index");
	check($res, "101112");
	
	//getValues
	$res = $db->getValues("SELECT TOP 5 Labornummer_Index FROM v_ngs_duo ORDER BY Labornummer_Index");
	check(count($res), 5);
}

end_test();

?>
