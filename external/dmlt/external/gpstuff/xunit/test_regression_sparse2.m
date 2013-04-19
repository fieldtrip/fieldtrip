function test_suite = test_regression_sparse2

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_REGRESSION_SPARSE2

initTestSuite;

function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_regression_sparse2')
demo_regression_sparse2
Eft_full = Eft_full(1:100);
Eft_var = Eft_var(1:100);
Varft_full = Varft_full(1:100);
Varft_var = Varft_var(1:100);
path = which('test_regression_sparse2.m');
path = strrep(path,'test_regression_sparse2.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testRegression_sparse2'); 
save(path, 'Eft_full', 'Eft_var', 'Varft_full', 'Varft_var');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictionsFull
values.real = load('realValuesRegression_sparse2.mat', 'Eft_full', 'Varft_full');
values.test = load(strrep(which('test_regression_sparse2.m'), 'test_regression_sparse2.m', 'testValues/testRegression_sparse2.mat'),'Eft_full', 'Varft_full');
assertElementsAlmostEqual(mean(values.real.Eft_full), mean(values.test.Eft_full), 'relative', 0.1);
assertElementsAlmostEqual(mean(values.real.Varft_full), mean(values.test.Varft_full), 'relative', 0.1);


function testPredictionsVar
values.real = load('realValuesRegression_sparse2.mat', 'Eft_var', 'Varft_var');
values.test = load(strrep(which('test_regression_sparse2.m'), 'test_regression_sparse2.m', 'testValues/testRegression_sparse2.mat'), 'Eft_var', 'Varft_var');
assertElementsAlmostEqual((values.real.Eft_var), (values.test.Eft_var), 'relative', 0.1);
assertElementsAlmostEqual((values.real.Varft_var), (values.test.Varft_var), 'relative', 0.1);
