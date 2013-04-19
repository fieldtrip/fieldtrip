function test_suite = test_regression_meanf

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_REGRESSION_MEANF

initTestSuite;



function testDemo
% Set random number stream so that test failing isn't because randomness.
% Run demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_regression_meanf')
demo_regression_meanf
path = which('test_regression_meanf.m');
path = strrep(path,'test_regression_meanf.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testRegression_meanf'); 
save(path, 'Eft');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all


function testPredictions
values.real = load('realValuesRegression_meanf.mat', 'Eft');
values.test = load(strrep(which('test_regression_meanf.m'), 'test_regression_meanf.m', 'testValues/testRegression_meanf.mat'), 'Eft');
assertElementsAlmostEqual(mean(values.real.Eft), mean(values.test.Eft), 'relative', 0.10);

