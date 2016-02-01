function test_suite = test_regression_additive2

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_REGRESSION_ADDITIVE2

initTestSuite;


function testDemo
   % Set random number stream so that test failing isn't because randomness.
   % Run demo & save test values.
   prevstream=setrandstream(0);
   
   disp('Running: demo_regression_additive2')
   demo_regression_additive2
    path = which('test_regression_additive2.m');
    path = strrep(path,'test_regression_additive2.m', 'testValues');
    if ~(exist(path, 'dir') == 7)
        mkdir(path)
    end
    path = strcat(path, '/testRegression_additive2');     
    save(path, 'Eft_map');
    
    % Set back initial random stream
    setrandstream(prevstream);
    drawnow;clear;close all
    
% Compare test values to real values.    

function testNeuralNetworkCFPrediction
    values.real = load('realValuesRegression_additive2.mat','Eft_map');
    values.test = load(strrep(which('test_regression_additive2.m'), 'test_regression_additive2.m', 'testValues/testRegression_additive2.mat'),'Eft_map');
    assertElementsAlmostEqual(mean(values.real.Eft_map), mean(values.test.Eft_map), 'relative', 0.05);