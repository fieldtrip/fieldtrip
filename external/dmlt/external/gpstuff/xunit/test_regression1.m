function test_suite = test_regression1

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_REGRESSION1

initTestSuite;


  function testDemo
    % Set random number stream so that test failing isn't because randomness.
    % Run demo & save test values.
    prevstream=setrandstream(0);
    
    disp('Running: demo_regression1')
    demo_regression1
    Eft_map = Eft_map(1:50);
    Varft_map = Varft_map(1:50);
    Eft_ia = Eft_ia(1:50);
    Varft_ia = Varft_ia(1:50);
    Eft_mc = Eft_mc(1:50);
    Varft_mc = Varft_mc(1:50);
    path = which('test_regression1.m');
    path = strrep(path,'test_regression1.m', 'testValues');
    if ~(exist(path, 'dir') == 7)
      mkdir(path)
    end
    path = strcat(path, '/testRegression1');
    save(path, 'K', 'C', 'w', 'Eft_map', 'Varft_map', ...
      'Eft_ia', 'Varft_ia', 'Eft_mc', 'Varft_mc');
    
    % Set back initial random stream
    setrandstream(prevstream);
    drawnow;clear;close all

% Test saved values with multiple tests. Covariance matrices and
% optimized parameters are tested with zero tolerance, while mean and 
% variances from various approximations are tested with relative tolerance
% (5% for MCMC and 1% for grid and IA)

    function testCovarianceMatrices
        values.real = load('realValuesRegression1.mat','K','C');
        values.test = load(strrep(which('test_regression1.m'), 'test_regression1.m', 'testValues/testRegression1.mat'),'K','C');
        assertElementsAlmostEqual(values.real.K, values.test.K);
        assertElementsAlmostEqual(values.real.C, values.test.C);
    

    function testOptimizedParameter
        values.real = load('realValuesRegression1.mat','w');
        values.test = load(strrep(which('test_regression1.m'), 'test_regression1.m', 'testValues/testRegression1.mat'),'w');
        assertElementsAlmostEqual(values.real.w, values.test.w, 'relative', 0.1);
        
    

    function testPredictedMeanVarianceGrid
        values.real = load('realValuesRegression1.mat','Eft_map','Varft_map');
        values.test = load(strrep(which('test_regression1.m'), 'test_regression1.m', 'testValues/testRegression1.mat'),'Eft_map','Varft_map');
        if length(values.test.Eft_map) > 50
            assertElementsAlmostEqual(values.real.Eft_map(1:50), values.test.Eft_map(1:50), 'relative', 0.1);
            assertElementsAlmostEqual(values.test.Varft_map(1:50), values.real.Varft_map(1:50), 'relative', 0.1);
        else
            assertElementsAlmostEqual(values.test.Eft_map, values.real.Eft_map, 'relative', 0.1);
            assertElementsAlmostEqual(values.test.Varft_map, values.real.Varft_map, 'relative', 0.1);
        end
        

    function testPredictedMeanVarianceMC
        values.real = load('realValuesRegression1.mat','Eft_mc','Varft_mc');
        values.test = load(strrep(which('test_regression1.m'), 'test_regression1.m', 'testValues/testRegression1.mat'),'Eft_mc','Varft_mc');
        if length(values.test.Eft_mc) > 50
            assertElementsAlmostEqual(mean(mean(values.test.Eft_mc(1:50))), mean(mean(values.real.Eft_mc(1:50))), 'relative', 0.1);
            assertElementsAlmostEqual(mean(values.test.Varft_mc(1:50)), mean(values.real.Varft_mc(1:50)), 'absolute', 0.3);
       else
            assertElementsAlmostEqual(mean(mean(values.test.Eft_mc)), mean(mean(values.real.Eft_mc)), 'relative', 0.1);
            assertElementsAlmostEqual(mean(mean(values.test.Varft_mc)), mean(mean(values.real.Varft_mc)), 'absolute', 0.3);
        end
   

    function testPredictedMeanVarianceIA
        values.real = load('realValuesRegression1.mat','Eft_ia','Varft_ia');
        values.test = load(strrep(which('test_regression1.m'), 'test_regression1.m', 'testValues/testRegression1.mat'),'Eft_ia','Varft_ia');
        if length(values.test.Eft_ia) > 50
            assertElementsAlmostEqual(mean(values.test.Eft_ia(1:50)), mean(values.real.Eft_ia(1:50)), 'relative', 0.1);
            assertElementsAlmostEqual(mean(values.test.Varft_ia(1:50)), mean(values.real.Varft_ia(1:50)), 'relative', 0.1);
        else
            assertElementsAlmostEqual(mean(values.test.Eft_ia), mean(values.real.Eft_ia), 'relative', 0.1);
            assertElementsAlmostEqual(mean(values.test.Varft_ia), mean(values.real.Varft_ia), 'relative', 0.1);
        end
    
    


