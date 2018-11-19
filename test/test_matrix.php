<?php

include("framework.php");

//##################################################################################
start_test("Matrix::slice");

$matrix = new Matrix();
$matrix->setHeaders(array("A", "B", "C", "D", "E", "F"));
$matrix->addRow(array(1 , 2, 3, 6, 5, 4));
$matrix->addRow(array(11,12,13,16,15,14));
$matrix->addRow(array(21,22,23,26,25,24));

$matrix->slice(array(1,3));

check($matrix->rows(), 3);
check($matrix->cols(), 2);
check($matrix->getHeaders(), array("B", "D"));
check($matrix->getRow(0), array(2, 6));
check($matrix->getRow(1), array(12, 16));
check($matrix->getRow(2), array(22, 26));

#part 2

$matrix = new Matrix();
$matrix->setHeaders(array("A", "B", "C", "D", "E", "F"));
$matrix->addRow(array(1 , 2, 3, 6, 5, 4));
$matrix->addRow(array(11,12,13,16,15,14));
$matrix->addRow(array(21,22,23,26,25,24));

$matrix->slice(array(2,3,4), true);

check($matrix->rows(), 3);
check($matrix->cols(), 6);
check($matrix->getHeaders(), array("C", "D", "E", "A", "B", "F"));
check($matrix->getRow(0), array( 3, 6, 5, 1, 2, 4));
check($matrix->getRow(1), array(13,16,15,11,12,14));
check($matrix->getRow(2), array(23,26,25,21,22,24));

end_test();

//##################################################################################
start_test("Matrix::cols");

$matrix = new Matrix();
check($matrix->cols(), 0);

end_test();

//##################################################################################
start_test("Matrix::rows");

$matrix = new Matrix();
check($matrix->rows(), 0);

end_test();

//##################################################################################
start_test("Matrix::comments");

$matrix = new Matrix();
check($matrix->comments(), 0);

end_test();

//##################################################################################
start_test("Matrix::addCol");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));

check($matrix->cols(), 1);
check($matrix->rows(), 6);

end_test();

//##################################################################################
start_test("Matrix::addRow");

$matrix = new Matrix();
$matrix->addRow(array(1,2,3,6,5,4));

check($matrix->cols(), 6);
check($matrix->rows(), 1);

end_test();

//##################################################################################
start_test("Matrix::get");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

check($matrix->get(1, 0), 2);
check($matrix->get(2, 1), "C");

end_test();


//##################################################################################
start_test("Matrix::getRow");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

check($matrix->getRow(1), array(2, "B"));

end_test();

//##################################################################################
start_test("Matrix::getCol");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

check($matrix->getCol(1), array("A", "B", "C", "D", "E", "F"));

end_test();

//##################################################################################
start_test("Matrix::getComments");

$matrix = new Matrix();

check($matrix->getComments(), array());

end_test();


//##################################################################################
start_test("Matrix::set");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

$matrix->set(2,1, "BLA");

check($matrix->get(2, 1), "BLA");

end_test();

//##################################################################################
start_test("Matrix::setCol");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

$matrix->setCol(0, array("X", "Y", "Z", "V", "W", "U"));
check($matrix->get(2, 0), "Z");
check($matrix->get(2, 1), "C");

end_test();

//##################################################################################
start_test("Matrix::setRow");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

$matrix->setRow(1, array("1A", "2B"));
check($matrix->get(2, 0), 3);
check($matrix->get(2, 1), "C");
check($matrix->get(1, 0), "1A");
check($matrix->get(1, 1), "2B");

end_test();


//##################################################################################
start_test("Matrix::setComments");

$matrix = new Matrix();
$matrix->setComments(array("a", "b"));
$comments = $matrix->getComments();

check(count($comments), 2);
check($comments[0], "a");
check($comments[1], "b");

end_test();

//##################################################################################
start_test("Matrix::addComment");

$matrix = new Matrix();
$matrix->addComment("a");
$comments = $matrix->getComments();
check(count($comments), 1);
check($comments[0], "a");

end_test();

//##################################################################################
start_test("Matrix::removeCol");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

$matrix->removeCol(0);

check($matrix->cols(), 1);
check($matrix->rows(), 6);
check($matrix->get(0, 0), "A");

end_test();

//##################################################################################
start_test("Matrix::removeRow");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

$matrix->removeRow(0);

check($matrix->cols(), 2);
check($matrix->rows(), 5);
check($matrix->get(0, 0), 2);
check($matrix->get(0, 1), "B");

end_test();


//##################################################################################
start_test("Matrix::getHeaders");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4), "a");
$matrix->addCol(array(1,2,3,6,5,4), "b");
$headers = $matrix->getHeaders();

check(count($headers), 2);
check($headers[0], "a");
check($headers[1], "b");

end_test();


//##################################################################################
start_test("Matrix::getHeader");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4), "a");
$matrix->addCol(array(1,2,3,6,5,4), "b");

check($matrix->getHeader(0), "a");
check($matrix->getHeader(1), "b");

end_test();

//##################################################################################
start_test("Matrix::setHeader");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4), "a");
$matrix->addCol(array(1,2,3,6,5,4), "b");
$matrix->setHeader(1, "c");
check($matrix->getHeader(0), "a");
check($matrix->getHeader(1), "c");

end_test();

//##################################################################################
start_test("Matrix::setHeaders");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4), "a");
$matrix->addCol(array(1,2,3,6,5,4), "b");
$matrix->setHeaders(array("c", "d"));
check($matrix->getHeader(0), "c");
check($matrix->getHeader(1), "d");

end_test();

//##################################################################################
start_test("Matrix::sort");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3,6,5,4));
$matrix->addCol(array("A", "B", "C", "D", "E", "F"));

$matrix->sort(0, SORT_NUMERIC);
check($matrix->getCol(0)==array(1,2,3,4,5,6), true);
check($matrix->getCol(1)==array("A", "B", "C", "F", "E", "D"), true);

$matrix->sort(1);
check($matrix->getCol(0)==array(1,2,3,6,5,4), true);
check($matrix->getCol(1)==array("A", "B", "C", "D", "E", "F"), true);

$matrix->sort(1, SORT_REGULAR, SORT_DESC);
check($matrix->getCol(0)==array(4,5,6,3,2,1), true);
check($matrix->getCol(1)==array("F", "E", "D", "C", "B", "A"), true);

end_test();

//##################################################################################
start_test("Matrix::fromTSV");

$matrix = Matrix::fromTSV(data_folder()."/matrix_fromTSV.txt");
check($matrix->cols(), 3);
check($matrix->rows(), 4);
check($matrix->comments(), 1);
check($matrix->getHeader(0), "a");
check($matrix->getHeader(1), "b");
check($matrix->getHeader(2), "c");
check($matrix->getCol(0), array(1, 4, 7, 10));
check($matrix->getCol(1), array(2, 5, 8, 11));
check($matrix->getCol(2), array(3, 6, 9, 12));
$comments = $matrix->getComments();
check($comments[0], "senseless comment");

$matrix = Matrix::fromTSV(data_folder()."/matrix_fromTSV.txt.gz");
check($matrix->cols(), 3);
check($matrix->rows(), 4);
check($matrix->comments(), 1);
check($matrix->getHeader(0), "a");
check($matrix->getHeader(1), "b");
check($matrix->getHeader(2), "c");
check($matrix->getCol(0), array(1, 4, 7, 10));
check($matrix->getCol(1), array(2, 5, 8, 11));
check($matrix->getCol(2), array(3, 6, 9, 12));
$comments = $matrix->getComments();
check($comments[0], "senseless comment");


$matrix = Matrix::fromTSV(data_folder()."/matrix_fromTSV_empty.txt");
check($matrix->cols(), 3);
check($matrix->rows(), 0);
check($matrix->comments(), 2);
check($matrix->getHeader(0), "a");
check($matrix->getHeader(1), "b");
check($matrix->getHeader(2), "c");

end_test();

//##################################################################################
start_test("Matrix::toTSV");

$filename = temp_file();
$matrix1 = Matrix::fromTSV(data_folder()."/matrix_toTSV.txt");
$matrix1->toTSV($filename);

$matrix = Matrix::fromTSV($filename);

check($matrix->cols(), 3);
check($matrix->rows(), 4);
check($matrix->comments(), 1);
check($matrix->getHeader(0), "a");
check($matrix->getHeader(1), "b");
check($matrix->getHeader(2), "c");
check($matrix->getCol(0), array(1, 4, 7, 10));
check($matrix->getCol(1), array(2, 5, 8, 11));
check($matrix->getCol(2), array(3, 6, 9, 12));
$comments = $matrix->getComments();
check($comments[0], "senseless comment");

end_test();

//##################################################################################
start_test("Matrix::transpose");

$matrix = new Matrix();
$matrix->addCol(array(1,2,3), "numbers");
$matrix->addCol(array("A", "B", "C"), "characters");

$matrix->transpose();
check($matrix->cols(), 3);
check($matrix->rows(), 2);
check($matrix->getCol(0), array(1, "A"));
check($matrix->getCol(1), array(2, "B"));
check($matrix->getCol(2), array(3, "C"));

$matrix->transpose();
check($matrix->cols(), 2);
check($matrix->rows(), 3);
check($matrix->getCol(0), array(1, 2, 3));
check($matrix->getCol(1), array("A", "B", "C"));

end_test();

//##################################################################################
start_test("Matrix::join");

$matrix = new Matrix();
$matrix->addCol(array(3,2,1), "numbers");
$matrix->addCol(array("C", "B", "A"), "characters");

$matrix2 = new Matrix();
$matrix2->addCol(array("A", "C", "B", "D"), "characters");
$matrix2->addCol(array(11, 33, 22, 44), "numbers");
$matrix2->addCol(array(111, 333, 222, 444), "numbers");

$matrix->join(1, $matrix2, 0);
check($matrix->cols(), 4);
check($matrix->getCol(0), array(3,2,1));
check($matrix->getCol(1), array("C", "B", "A"));
check($matrix->getCol(2), array(33,22,11));
check($matrix->getCol(3), array(333,222,111));

//##################################################################################
start_test("Matrix::join");

$matrix = new Matrix();
$matrix->addCol(array(3,2,1), "numbers");
$matrix->addCol(array("C", "B", "A"), "characters");

$matrix2 = new Matrix();
$matrix2->addCol(array("A", "C", "B", "D"), "characters");
$matrix2->addCol(array(11, 33, 22, 44), "numbers");
$matrix2->addCol(array(111, 333, 222, 444), "numbers");

$matrix->join(1, $matrix2, 0);
check($matrix->cols(), 4);
check($matrix->getCol(0), array(3,2,1));
check($matrix->getCol(1), array("C", "B", "A"));
check($matrix->getCol(2), array(33,22,11));
check($matrix->getCol(3), array(333,222,111));

end_test();


//##################################################################################
//still need to improve matrix-test
start_test("Matrix::resize");

$matrix = new Matrix();
$matrix->resize(10, 10);
check($matrix->rows(), 10);
check($matrix->cols(), 10);
$matrix->resize(10, 8);
check($matrix->rows(), 10);
check($matrix->cols(), 8);
$matrix->resize(8, 8);
check($matrix->rows(), 8);
check($matrix->cols(), 8);
$matrix->resize(6, 6);
check($matrix->rows(), 6);
check($matrix->cols(), 6);
$matrix->resize(9, 3);
check($matrix->rows(), 9);
check($matrix->cols(), 3);
$matrix->resize(2, 9);
check($matrix->rows(), 2);
check($matrix->cols(), 9);
$matrix->resize(0, 0);
check($matrix->rows(), 0);
check($matrix->cols(), 0);

end_test();

//##################################################################################
start_test("Matrix::unique");

$matrix1 = new Matrix();
$matrix1->addRow(array("1","2","3","4"));
$matrix1->addRow(array("1","3","3","4"));
$matrix1->addRow(array("1","2","3","4"));
$matrix1->addRow(array("1","4","3","4"));
$matrix1->unique();

$matrix2 = new Matrix();
$matrix2->addRow(array("1","2","3","4"));
$matrix2->addRow(array("1","3","3","4"));
$matrix2->addRow(array("1","4","3","4"));

check(serialize($matrix1), serialize($matrix2));

end_test();


//##################################################################################
start_test("Matrix::getColumnIndex");

$matrix = new Matrix();
$matrix->setHeaders(array("c1","c2","c3","c4"));
$matrix->addRow(array("1","2","3","4"));
$matrix->addRow(array("1","3","3","4"));
$matrix->addRow(array("3","2","3","4"));
$matrix->addRow(array("1","4","3","4"));

check($matrix->getColumnIndex("c3"), 2);
check($matrix->getColumnIndex("c1"), 0);

check($matrix->getColumnIndex("c4", true), 3);
check($matrix->getColumnIndex("c2", true), 1);

check($matrix->getColumnIndex("c0", false, false), false);
check($matrix->getColumnIndex("c5", false, false), false);

end_test();

//##################################################################################
start_test("Matrix::getData");

$matrix = new Matrix();
$matrix->addRow(array(1,2,3,6,5,4));
$matrix->addRow(array("A", "B", "C", "D", "E", "F"));

check($matrix->getData(), array(array(1,2,3,6,5,4), array("A", "B", "C", "D", "E", "F")));


end_test();

//##################################################################################
start_test("Matrix::insertCol");

$matrix = new Matrix();
$matrix->setHeaders(array("c1","c2","c3","c4"));
$matrix->addRow(array("1","2","3","4"));
$matrix->addRow(array("1","3","3","4"));
$matrix->addRow(array("3","2","3","4"));
$matrix->addRow(array("1","4","3","4"));

$matrix->insertCol(1, array("5","5","5","5"));

check($matrix->get(1,0), "1");
check($matrix->get(1,1), "5");
check($matrix->get(1,2), "3");
check($matrix->get(3,0), "1");
check($matrix->get(3,1), "5");
check($matrix->get(3,2), "4");

end_test();
?>
