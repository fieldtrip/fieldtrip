function test_suite = test_spatial1

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_SPATIAL1

initTestSuite;


function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_spatial1')
demo_spatial1
Ef = Ef(1:100);
Varf = Varf(1:100);
path = which('test_spatial1.m');
path = strrep(path,'test_spatial1.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testSpatial1'); 
save(path, 'Elth', 'Elth2', 'Ef', 'Varf');

% Set back initial random stream
setrandstream(prevstream);

% Compare test values to real values.

function testEstimatesIA
values.real = load('realValuesSpatial1.mat', 'Elth', 'Elth2');
values.test = load(strrep(which('test_spatial1.m'), 'test_spatial1.m', 'testValues/testSpatial1.mat'), 'Elth', 'Elth2');
assertElementsAlmostEqual(values.real.Elth, values.test.Elth, 'relative', 0.1);
assertElementsAlmostEqual(values.real.Elth2, values.test.Elth2, 'relative', 0.1);


function testPredictionIA
values.real = load('realValuesSpatial1.mat', 'Ef', 'Varf');
values.test = load(strrep(which('test_spatial1.m'), 'test_spatial1.m', 'testValues/testSpatial1.mat'), 'Ef', 'Varf');
assertElementsAlmostEqual(mean(values.real.Ef), mean(values.test.Ef), 'relative', 0.1);
assertElementsAlmostEqual(mean(values.real.Varf), mean(values.test.Varf), 'relative', 0.1);

