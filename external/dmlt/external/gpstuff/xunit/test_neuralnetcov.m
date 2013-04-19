function test_suite = test_neuralnetcov

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_NEURALNETCOV

initTestSuite;


function testDemo
% Set random number stream so that failing isn't because randomness. Run
% demo & save test values.
prevstream=setrandstream(0);

disp('Running: demo_neuralnetcov')
demo_neuralnetcov
path = which('test_neuralnetcov.m');
path = strrep(path,'test_neuralnetcov.m', 'testValues');
if ~(exist(path, 'dir') == 7)
    mkdir(path)
end
path = strcat(path, '/testNeuralnetcov'); 
save(path,  'Eft_map', 'Varft_map', 'Eft_map2', 'Varft_map2');

% Set back initial random stream
setrandstream(prevstream);
drawnow;clear;close all

% Compare test values to real values.

function testPredictions
values.real = load('realValuesNeuralnetcov.mat','Eft_map', 'Eft_map2','Varft_map','Varft_map2');
values.test = load(strrep(which('test_neuralnetcov.m'), 'test_neuralnetcov.m', 'testValues/testNeuralnetcov.mat'), 'Eft_map', 'Eft_map2','Varft_map','Varft_map2');
assertElementsAlmostEqual(mean(values.real.Eft_map), mean(values.test.Eft_map), 'relative', 0.05);
assertElementsAlmostEqual(mean(values.real.Eft_map2), mean(values.test.Eft_map2), 'relative', 0.05);
assertElementsAlmostEqual(mean(values.real.Varft_map), mean(values.test.Varft_map), 'relative', 0.05);
assertElementsAlmostEqual(mean(values.real.Varft_map2), mean(values.test.Varft_map2), 'relative', 0.05);
