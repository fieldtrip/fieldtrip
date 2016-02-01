function test_suite = test_classific

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
    
    disp('Running: demo_classific')
    demo_classific
    Eft_la = Eft_la(1:100);
    Eft_ep = Eft_ep(1:100);
    Efs_mc = Efs_mc(1:50,1:50);
    Varft_la = Varft_la(1:100);
    Varft_ep = Varft_ep(1:100);
    Varfs_mc = Varfs_mc(1:50,1:50);
    path = which('test_classific.m');
    path = strrep(path,'test_classific.m', 'testValues');
    if ~(exist(path, 'dir') == 7)
      mkdir(path)
    end
    path = strcat(path, '/testClassific');
    save(path, 'Eft_la', 'Varft_la', 'Eft_ep', 'Varft_ep', 'Efs_mc', 'Varfs_mc');
    
    % Set back initial random stream
    setrandstream(prevstream);
    drawnow;clear;close all

 % Compare test values to real values.
        
  function testLaplace
    values.real = load('realValuesClassific.mat','Eft_la','Varft_la');
    values.test = load(strrep(which('test_classific.m'), 'test_classific.m', 'testValues/testClassific.mat'),'Eft_la','Varft_la');
    assertElementsAlmostEqual(mean(values.real.Eft_la),mean(values.test.Eft_la),'relative',0.01);
    assertElementsAlmostEqual(mean(values.real.Varft_la),mean(values.test.Varft_la),'relative',0.01);
        
        
  function testEP
    values.real = load('realValuesClassific.mat','Eft_ep','Varft_ep');
    values.test = load(strrep(which('test_classific.m'), 'test_classific.m', 'testValues/testClassific.mat'),'Eft_ep','Varft_ep');
    assertElementsAlmostEqual(mean(values.real.Eft_ep),mean(values.test.Eft_ep),'relative',0.01);
    assertElementsAlmostEqual(mean(values.real.Varft_ep),mean(values.test.Varft_ep),'relative',0.01);


  function testMC
    values.real = load('realValuesClassific.mat','Efs_mc','Varfs_mc');
    values.test = load(strrep(which('test_classific.m'), 'test_classific.m', 'testValues/testClassific.mat'),'Efs_mc','Varfs_mc');
    assertElementsAlmostEqual(mean(mean(values.real.Efs_mc)),mean(mean(values.test.Efs_mc)),'absolute',0.7);
    assertElementsAlmostEqual(mean(mean(values.real.Varfs_mc)),mean(mean(values.test.Varfs_mc)),'absolute',0.7);