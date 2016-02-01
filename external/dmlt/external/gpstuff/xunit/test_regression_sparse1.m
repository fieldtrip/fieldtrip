function test_suite = test_regression_sparse1

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_REGRESSION_SPARSE1

initTestSuite;


function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_regression_sparse1')
demo_regression_sparse1
path = which('test_regression_sparse1');
path = strrep(path,'test_regression_sparse1.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testRegression_sparse1'); 
save(path, 'Eft_fic', 'Eft_pic', 'Eft_var', 'Eft_dtc', 'Eft_cs');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictionsCS
values.real = load('realValuesRegression_sparse1', 'Eft_cs');
values.test = load(strrep(which('test_regression_sparse1.m'), 'test_regression_sparse1.m', 'testValues/testRegression_sparse1'), 'Eft_cs');
assertElementsAlmostEqual((values.real.Eft_cs), (values.test.Eft_cs), 'absolute', 0.1);


function testPredictionsFIC
values.real = load('realValuesRegression_sparse1', 'Eft_fic');
values.test = load(strrep(which('test_regression_sparse1.m'), 'test_regression_sparse1.m', 'testValues/testRegression_sparse1'), 'Eft_fic');
assertElementsAlmostEqual((values.real.Eft_fic), (values.test.Eft_fic), 'absolute', 0.1);


function testPredictionsPIC
values.real = load('realValuesRegression_sparse1', 'Eft_pic');
values.test = load(strrep(which('test_regression_sparse1.m'), 'test_regression_sparse1.m', 'testValues/testRegression_sparse1'), 'Eft_pic');
assertElementsAlmostEqual((values.real.Eft_pic), (values.test.Eft_pic), 'absolute', 0.1);


function testPredictionsVAR
values.real = load('realValuesRegression_sparse1', 'Eft_var');
values.test = load(strrep(which('test_regression_sparse1.m'), 'test_regression_sparse1.m', 'testValues/testRegression_sparse1'), 'Eft_var');
assertElementsAlmostEqual((values.real.Eft_var), (values.test.Eft_var), 'absolute', 0.1);


function testPredictionsDTC
values.real = load('realValuesRegression_sparse1', 'Eft_dtc');
values.test = load(strrep(which('test_regression_sparse1.m'), 'test_regression_sparse1.m', 'testValues/testRegression_sparse1'), 'Eft_dtc');
assertElementsAlmostEqual((values.real.Eft_dtc), (values.test.Eft_dtc), 'absolute', 0.1);
