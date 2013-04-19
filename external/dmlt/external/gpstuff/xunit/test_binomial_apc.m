function test_suite = test_binomial_apc

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_BINOMIAL_APC

% Copyright (c) 2011-2012 Ville Tolvanen

initTestSuite;


  function testDemo
    % Set random number stream so that the test failing isn't because
    % randomness. Run demo & save test values
    prevstream=setrandstream(0);
    
    disp('Running: demo_binomial_apc')
    demo_binomial_apc
    path = which('test_binomial_apc.m');
    path = strrep(path,'test_binomial_apc.m', 'testValues');
    if ~(exist(path, 'dir') == 7)
      mkdir(path)
    end
    path = strcat(path, '/testBinomial_apc');
    save(path, 'Eft', 'Varft', 'Eft_3', 'Varft_3');
    
    % Set back initial random stream
    setrandstream(prevstream);
    drawnow;clear;close all
    
    % Compare test values to real values.
    
  function testPredictionsAll
    values.real = load('realValuesBinomial_apc', 'Eft', 'Varft');
    values.test = load(strrep(which('test_binomial_apc.m'), 'test_binomial_apc.m', 'testValues/testBinomial_apc'), 'Eft', 'Varft');
    assertElementsAlmostEqual((values.real.Eft), (values.test.Eft), 'absolute', 0.1);
    assertElementsAlmostEqual((values.real.Varft), (values.test.Varft), 'absolute', 0.1);
    
  function testPredictionsCohort
    values.real = load('realValuesBinomial_apc', 'Eft_3', 'Varft_3');
    values.test = load(strrep(which('test_binomial_apc.m'), 'test_binomial_apc.m', 'testValues/testBinomial_apc'), 'Eft_3', 'Varft_3');
    assertElementsAlmostEqual((values.real.Eft_3), (values.test.Eft_3), 'absolute', 0.1);
    assertElementsAlmostEqual((values.real.Varft_3), (values.test.Varft_3), 'absolute', 0.1);
