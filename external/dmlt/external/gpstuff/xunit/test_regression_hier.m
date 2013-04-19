function test_suite = test_regression_hier

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_REGRESSION_HIER

initTestSuite;

% Set random number stream so that test failing isn't because randomness.
% Run demo & save test values.

function testDemo
prevstream=setrandstream(0);

disp('Running: demo_regression_hier')
demo_regression_hier
path = which('test_regression_hier.m');
path = strrep(path,'test_regression_hier.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testRegression_hier'); 
save(path, 'Eff');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all


function testPredictionMissingData
values.real = load('realValuesRegression_hier', 'Eff');
values.test = load(strrep(which('test_regression_hier.m'), 'test_regression_hier.m', 'testValues/testRegression_hier'), 'Eff');
assertVectorsAlmostEqual(mean(values.real.Eff), mean(values.test.Eff), 'relative', 0.01);

