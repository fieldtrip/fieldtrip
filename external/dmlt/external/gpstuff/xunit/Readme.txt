Tests in this directory use MATLAB xUnit Test Framework toolbox.
Test structure in every test is basically the same: Run demo, save values from demo, compare those saved test values to previously saved "correct" values. Users can run these tests to determine if some calculations within demos don't give correct answers.
Every m-file in the directory is written so that it includes(or can include) multiple tests for a given demo. There is one m-file for every demo made while writing this document. 
Script test_all compiles all the tests in the directory in one TestSuite object and runs them. After running the tests, logger object gives test results, what went wrong etc.
As previously mentioned, every test in the directory, present at the moment, uses same build for testing. First, to include multiple tests in one m-file, you should start your own test code with these two lines:

1	function test_suite = testName
2	initTestSuite;

Every test file and test should include string "test", either in front of the name or after the name, e.g. testName or nameTest. To write the actual tests in the m-file, you use subfunctions, including string "test" again.

4	function testDemo
5		demo_regression1;
6		save('testValues/testRegression1');
7		drawnow;clear;close all

xUnit doesn't recognize tests if the name of the main m-file, or the subfunction, doesn't have string "test" in it. For transparency, every test written at the moment, is named test_demoName. 
To compare the test values to real values, most of the tests use two test functions assertElementsAlmostEqual and assertVectorsAlmostEqual, provided by the xUnit Test Framework. 
Function assertElementsAlmostEqual(a, b, options) compares elements of a to b, with 'options' providing method and tolerance level, e.g. assertElementsAlmostEqual(1, 2, 'absolute', 0.5) returns error, 
while assertElementsAlmostEqual(1, 1.2, 'relative', 0.2) doesn't. By default, the comparison uses relative tolerance test. Function assertVectorsAlmostEqual works the same way. 
While writing your own test, it might also be wise to set the random number stream before taking real values and test values so that the test failing isn't because randomness.

9	function testElements
10		k = 10.1
11		l = 10
12		assertElementsAlmostEqual(k, l, 'absolute', 0.5)  % Works fine
13		assertElementsAlmostEqual(k, l, 'absolute', 0.001)  % Error
14		assertElementsAlmostEqual(k, l, 'relative', 0.01) % Works fine(1% error tolerance level)


You can run a specific test with command "runtests testName", and all the tests in the root directory(functions which include string "test") simply by writing "runtests". 
As mentioned before, script testAll runs all of the tests also, but with it, you can check results later, simply by examining logger object.

For more thorough tutorial to MATLAB xUnit Test Framework and for some additional testing functions, visit 
http://www.mathworks.com/matlabcentral/fx_files/22846/11/content/matlab_xunit/doc/xunit_product_page.html

Real values from revision 990.


