function test_suite = test_loopred

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_BINOMIAL1

% Copyright (c) 2011-2012 Ville Tolvanen

initTestSuite;


  function testDemo
    % Set random number stream so that the test failing isn't because
    % randomness. Run demo & save test values.
    prevstream=setrandstream(0);    
    disp('Running: demo_loopred')
    demo_loopred
    path = which('test_loopred');
    path = strrep(path,'test_loopred.m', 'testValues');
    if ~(exist(path, 'dir') == 7)
      mkdir(path)
    end
    path = strcat(path, '/testLoopred');
    save(path, 'Eft_lrs', 'Varft_lrs', 'lpyt_lrs', 'Eft_cav', 'Varft_cav', 'lpyt_cav', 'Eft_inla', 'Varft_inla', 'lpyt_inla' ...
      , 'Eft_ep', 'Varft_ep', 'lpyt_ep');
    
    % Set back initial random stream
    setrandstream(prevstream);
    drawnow;clear;close all


% Test predictive mean, variance and density for binomial model with 5% tolerance.        
        
  function testLrsLOO
    values.real = load('realValuesLoopred.mat','Eft_lrs','Varft_lrs','lpyt_lrs');
    values.test = load(strrep(which('test_loopred.m'), 'test_loopred.m', 'testValues/testLoopred.mat'),'Eft_lrs','Varft_lrs','lpyt_lrs');
    assertElementsAlmostEqual(values.real.Eft_lrs, values.test.Eft_lrs, 'absolute', 0.1);
    assertElementsAlmostEqual(values.real.Varft_lrs, values.test.Varft_lrs, 'absolute', 0.1);
    assertElementsAlmostEqual(values.real.lpyt_lrs, values.test.lpyt_lrs, 'absolute', 0.1);

  function testCavityLOO
    values.real = load('realValuesLoopred.mat','Eft_cav','Varft_cav','lpyt_cav');
    values.test = load(strrep(which('test_loopred.m'), 'test_loopred.m', 'testValues/testLoopred.mat'),'Eft_cav','Varft_cav','lpyt_cav');
    assertElementsAlmostEqual(values.real.Eft_cav, values.test.Eft_cav, 'absolute', 0.1);
    assertElementsAlmostEqual(values.real.Varft_cav, values.test.Varft_cav, 'absolute', 0.1);
    assertElementsAlmostEqual(values.real.lpyt_cav, values.test.lpyt_cav, 'absolute', 0.1);
    
    
  function testInlaLOO
    values.real = load('realValuesLoopred.mat','Eft_inla','Varft_inla','lpyt_inla');
    values.test = load(strrep(which('test_loopred.m'), 'test_loopred.m', 'testValues/testLoopred.mat'),'Eft_inla','Varft_inla','lpyt_inla');
    assertElementsAlmostEqual(values.real.Eft_inla, values.test.Eft_inla, 'absolute', 0.1);
    assertElementsAlmostEqual(values.real.Varft_inla, values.test.Varft_inla, 'absolute', 0.1);
    assertElementsAlmostEqual(values.real.lpyt_inla, values.test.lpyt_inla, 'absolute', 0.1); 
    
  function testEPLOO
    values.real = load('realValuesLoopred.mat','Eft_ep','Varft_ep','lpyt_ep');
    values.test = load(strrep(which('test_loopred.m'), 'test_loopred.m', 'testValues/testLoopred.mat'),'Eft_ep','Varft_ep','lpyt_ep');
    assertElementsAlmostEqual(values.real.Eft_ep, values.test.Eft_ep, 'absolute', 0.1);
    assertElementsAlmostEqual(values.real.Varft_ep, values.test.Varft_ep, 'absolute', 0.1);
    assertElementsAlmostEqual(values.real.lpyt_ep, values.test.lpyt_ep, 'absolute', 0.1);    
 

