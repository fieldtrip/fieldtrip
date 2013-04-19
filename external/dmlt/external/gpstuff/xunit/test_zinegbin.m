function test_suite = test_zinegbin

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_ZINEGBIN

% Copyright (c) 2011-2012 Ville Tolvanen

initTestSuite;

function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_zinegbin')
demo_zinegbin;
path = which('test_zinegbin.m');
path = strrep(path,'test_zinegbin.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testZinegbin'); 
Ef=Ef(1:100); Varf=diag(Varf(1:100,1:100));
save(path, 'Ef', 'Varf');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictionsZinegbin
values.real = load('realValuesZinegbin', 'Ef', 'Varf');
values.test = load(strrep(which('test_zinegbin.m'), 'test_zinegbin.m', 'testValues/testZinegbin'), 'Ef', 'Varf');
assertElementsAlmostEqual(values.real.Ef, values.test.Ef, 'relative', 0.10);
assertElementsAlmostEqual(values.real.Varf, values.test.Varf, 'relative', 0.10);
