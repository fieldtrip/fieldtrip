function test_suite = test_spatial2

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_SPATIAL2

initTestSuite;



function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_spatial2')
demo_spatial2
Ef = Ef(1:100);
Varf = Varf(1:100);
C = C(1:50, 1:50);
path = which('test_spatial2.m');
path = strrep(path,'test_spatial2.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testSpatial2'); 
save(path, 'Ef', 'Varf', 'C');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictionsEP
values.real = load('realValuesSpatial2.mat', 'Ef', 'Varf');
values.test = load(strrep(which('test_spatial2.m'), 'test_spatial2.m', 'testValues/testSpatial2.mat'), 'Ef', 'Varf');
assertElementsAlmostEqual(mean(values.test.Ef), mean(values.real.Ef), 'relative', 0.1);
assertElementsAlmostEqual(mean(values.test.Varf), mean(values.real.Varf), 'relative', 0.1);


function testCovarianceMatrix
values.real = load('realValuesSpatial2.mat', 'C');
values.test = load(strrep(which('test_spatial2.m'), 'test_spatial2.m', 'testValues/testSpatial2.mat'), 'C');
assertElementsAlmostEqual(mean(values.real.C), mean(values.test.C), 'relative', 0.1);