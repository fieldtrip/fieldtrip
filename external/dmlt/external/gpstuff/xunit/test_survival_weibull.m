function test_suite = test_survival_weibull

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_SURVIVAL_WEIBULL

% Copyright (c) 2011-2012 Ville Tolvanen

initTestSuite;

function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_survival_weibull')
demo_survival_weibull;
path = which('test_survival_weibull.m');
path = strrep(path,'test_survival_weibull.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testSurvival_weibull'); 
save(path, 'Ef1', 'Ef2', 'Varf1', 'Varf2');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictionsWeibull
values.real = load('realValuesSurvival_weibull', 'Ef1', 'Varf1', 'Ef2', 'Varf2');
values.test = load(strrep(which('test_survival_weibull.m'), 'test_survival_weibull.m', 'testValues/testSurvival_weibull'), 'Ef1', 'Varf1', 'Ef2', 'Varf2');
assertElementsAlmostEqual(values.real.Ef1, values.test.Ef1, 'relative', 0.10);
assertElementsAlmostEqual(values.real.Ef2, values.test.Ef2, 'relative', 0.10);
assertElementsAlmostEqual(values.real.Varf1, values.test.Varf1, 'relative', 0.10);
assertElementsAlmostEqual(values.real.Varf2, values.test.Varf2, 'relative', 0.10);