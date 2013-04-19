function test_suite = test_survival_coxph

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_SURVIVAL_COXPH

% Copyright (c) 2011-2012 Ville Tolvanen

initTestSuite;

function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_survival_coxph')
demo_survival_coxph;
path = which('test_survival_coxph.m');
path = strrep(path,'test_survival_coxph.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testSurvival_coxph'); 
save(path, 'Ef1', 'Ef2', 'Varf1', 'Varf2');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictionsCoxph
values.real = load('realValuesSurvival_coxph', 'Ef1', 'Varf1', 'Ef2', 'Varf2');
values.test = load(strrep(which('test_survival_coxph.m'), 'test_survival_coxph.m', 'testValues/testSurvival_coxph'), 'Ef1', 'Varf1', 'Ef2', 'Varf2');
assertElementsAlmostEqual(values.real.Ef1, values.test.Ef1, 'absolute', 0.10);
assertElementsAlmostEqual(values.real.Ef2, values.test.Ef2, 'absolute', 0.10);
assertElementsAlmostEqual(values.real.Varf1, values.test.Varf1, 'absolute', 0.10);
assertElementsAlmostEqual(values.real.Varf2, values.test.Varf2, 'absolute', 0.10);