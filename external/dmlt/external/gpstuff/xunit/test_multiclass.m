function test_suite = test_multiclass

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_MULTICLASS

initTestSuite;


function testDemo
    % Set random number stream so that test failing isn't because randomness.
    % Run demo & save test values.
    prevstream=setrandstream(0);
    
    disp('Running: demo_multiclass')
    demo_multiclass
    Eft=Eft(1:100,1:3);
    Varft=Varft(1:3,1:3,1:100);
    Covft=Covft(1:3,1:3,1:100);
    path = which('test_multiclass.m');
    path = strrep(path,'test_multiclass.m', 'testValues');
    if ~(exist(path, 'dir') == 7)
        mkdir(path)
    end
    path = strcat(path, '/testMulticlass');     
    save(path,'Eft','Varft','Covft');
    
    % Set back initial random stream
    setrandstream(prevstream);
    drawnow;clear;close all
    
% Compare test values to real values.

function testPredictions
    values.real = load('realValuesMulticlass.mat');
    values.test = load(strrep(which('test_multiclass.m'), 'test_multiclass.m', 'testValues/testMulticlass.mat'));
    assertElementsAlmostEqual(mean(values.real.Eft), mean(values.test.Eft), 'relative', 0.01);
    assertElementsAlmostEqual(mean(values.real.Varft), mean(values.test.Varft), 'relative', 0.01);
    assertElementsAlmostEqual(mean(values.real.Covft), mean(values.test.Covft), 'relative', 0.01);
