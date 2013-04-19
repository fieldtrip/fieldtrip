function test_suite = test_regression_robust

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_REGRESSION_ROBUST

initTestSuite;


function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_regression_robust')
demo_regression_robust
path = which('test_regression_robust.m');
path = strrep(path,'test_regression_robust.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testRegression_robust'); 
w=gp_pak(rr);
save(path, 'Eft', 'Varft', 'w')

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictionEP
values.real = load('realValuesRegression_robust', 'Eft', 'Varft');
values.test = load(strrep(which('test_regression_robust.m'), 'test_regression_robust.m', 'testValues/testRegression_robust'), 'Eft', 'Varft');
assertElementsAlmostEqual(mean(values.real.Eft), mean(values.test.Eft), 'relative', 0.05);
assertElementsAlmostEqual(mean(values.real.Varft), mean(values.test.Varft), 'relative', 0.05);

function testMCMCSamples
values.real = load('realValuesRegression_robust', 'w');
values.test = load(strrep(which('test_regression_robust.m'), 'test_regression_robust.m', 'testValues/testRegression_robust'), 'w');
assertElementsAlmostEqual(mean(values.real.w), mean(values.test.w), 'relative', 0.01);
assertElementsAlmostEqual(mean(values.real.w), mean(values.test.w), 'relative', 0.01);

