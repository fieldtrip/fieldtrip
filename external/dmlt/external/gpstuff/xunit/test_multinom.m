function test_suite = test_multinom

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

disp('Running: demo_multinom')
demo_multinom;
path = which('test_multinom.m');
path = strrep(path,'test_multinom.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testMultinom'); 
save(path, 'Eft', 'pyt2');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictionsMultinom
values.real = load('realValuesMultinom', 'Eft', 'pyt2');
values.test = load(strrep(which('test_multinom.m'), 'test_multinom.m', 'testValues/testMultinom'), 'Eft', 'pyt2');
assertElementsAlmostEqual(values.real.Eft, values.test.Eft, 'absolute', 0.10);
assertElementsAlmostEqual(values.real.pyt2, values.test.pyt2, 'absolute', 0.10);

