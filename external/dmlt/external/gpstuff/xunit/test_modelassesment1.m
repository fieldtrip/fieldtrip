function test_suite = test_modelassesment1

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_SURVIVAL_WEIBULL

% Copyright (c) 2011-2012 Ville Tolvanen

initTestSuite;


  function testDemo
    % Set random number stream so that test failing isn't because randomness.
    % Run demo % save test values.
    prevstream=setrandstream(0);
    
    disp('Running: demo_modelassesment1')
    demo_modelassesment1
    path = which('test_modelassesment1.m');
    path = strrep(path,'test_modelassesment1.m', 'testValues');
    if ~(exist(path, 'dir') == 7)
      mkdir(path)
    end
    path = strcat(path, '/testModelAssesment1');
    save(path, 'DIC', 'DIC2', 'DIC_latent', ...
      'p_eff', 'p_eff2', 'p_eff_latent', 'p_eff_latent2', 'mlpd_cv', 'rmse_cv', 'mlpd_loo', 'rmse_loo');
    
    % Set back initial random stream
    setrandstream(prevstream);
    drawnow;clear;close all
        
  % Compare test values to real values.

  function testDICParameters
    values.real = load('realValuesModelAssesment1.mat','DIC'); 
    values.real.DIC(isnan(values.real.DIC)) = 0;
    values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'),'DIC');
    values.test.DIC(isnan(values.test.DIC)) = 0;
    assertVectorsAlmostEqual(values.real.DIC, values.test.DIC, 'relative', 0.20);
    
        
  function testDICAll
    values.real = load('realValuesModelAssesment1.mat','DIC2');
    values.real.DIC2(isnan(values.real.DIC2)) = 0;
    values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'),'DIC2');
    values.test.DIC2(isnan(values.test.DIC2)) = 0;
    assertVectorsAlmostEqual(values.real.DIC2, values.test.DIC2, 'relative', 0.20);
    
        
  function testDICLatent
    values.real = load('realValuesModelAssesment1.mat','DIC_latent');
    values.real.DIC_latent(isnan(values.real.DIC_latent)) = 0;
    values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'),'DIC_latent');
    values.test.DIC_latent(isnan(values.test.DIC_latent)) = 0;
    assertVectorsAlmostEqual(values.real.DIC_latent, values.test.DIC_latent, 'relative', 0.20);
    
        
   function testPeffLatentMarginalized
     values.real = load('realValuesModelAssesment1.mat','p_eff');
     values.real.p_eff(isnan(values.real.p_eff)) = 0;
     values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'),'p_eff');
     values.test.p_eff(isnan(values.test.p_eff)) = 0;
     assertVectorsAlmostEqual(values.real.p_eff, values.test.p_eff, 'relative', 0.20);
     

   function testPeffAll
     values.real = load('realValuesModelAssesment1.mat','p_eff2');
     values.real.p_eff2(isnan(values.real.p_eff2)) = 0;
     values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'),'p_eff2');
     values.test.p_eff2(isnan(values.test.p_eff2)) = 0;
     assertVectorsAlmostEqual(values.real.p_eff2, values.test.p_eff2, 'relative', 0.20);
     
        
   function testPeffLatent
     values.real = load('realValuesModelAssesment1.mat','p_eff_latent');
     values.real.p_eff_latent(isnan(values.real.p_eff_latent)) = 0;
     values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'),'p_eff_latent');
     values.test.p_eff_latent(isnan(values.test.p_eff_latent)) = 0;
     assertVectorsAlmostEqual(values.real.p_eff_latent, values.test.p_eff_latent, 'relative', 0.20);
     
        
   function testPeffLatent2
     values.real = load('realValuesModelAssesment1.mat','p_eff_latent2');
     values.real.p_eff_latent2(isnan(values.real.p_eff_latent2)) = 0;
     values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'),'p_eff_latent2');
     values.test.p_eff_latent2(isnan(values.test.p_eff_latent2)) = 0;
     assertVectorsAlmostEqual(values.real.p_eff_latent2, values.test.p_eff_latent2, 'relative', 0.20);
     
     
   function testLogPredDensity10foldCV
     values.real = load('realValuesModelAssesment1.mat', 'mlpd_cv');
     values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'), 'mlpd_cv');
     assertVectorsAlmostEqual(values.real.mlpd_cv, values.test.mlpd_cv, 'relative', 0.20);
     
        
   function testMeanSquaredError10foldCV
     values.real = load('realValuesModelAssesment1.mat', 'rmse_cv');
     values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'), 'rmse_cv');
     assertVectorsAlmostEqual(values.real.rmse_cv, values.test.rmse_cv, 'relative', 0.20);
                 
   function testLogPredDensityLOOPRED
     values.real = load('realValuesModelAssesment1.mat', 'mlpd_loo');
     values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'), 'mlpd_loo');
     assertVectorsAlmostEqual(values.real.mlpd_loo, values.test.mlpd_loo, 'relative', 0.20);
     
        
   function testMeanSquaredErrorLOOPRED
     values.real = load('realValuesModelAssesment1.mat', 'rmse_loo');
     values.test = load(strrep(which('test_modelassesment1.m'), 'test_modelassesment1.m', 'testValues/testModelAssesment1.mat'), 'rmse_loo');
     assertVectorsAlmostEqual(values.real.rmse_loo, values.test.rmse_loo, 'relative', 0.20);
                 