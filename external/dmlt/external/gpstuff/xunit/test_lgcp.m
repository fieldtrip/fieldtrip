function test_suite = test_lgcp

%   Run specific demo and save values for comparison.
%
%   See also
%     TEST_ALL, DEMO_LGCP

% Copyright (c) 2011-2012 Ville Tolvanen

initTestSuite;


  function testDemo
    % Set random number stream so that test failing isn't because randomness.
    % Run demo & save test values.
    prevstream=setrandstream(0);
    
    disp('Running: demo_lgcp')
    demo_lgcp
    path = which('test_lgcp.m');
    path = strrep(path,'test_lgcp.m', 'testValues');
    if ~(exist(path, 'dir') == 7)
      mkdir(path)
    end
    path = strcat(path, '/testLgcp');
    save(path)
    
    % Set back initial random stream
    setrandstream(prevstream);
    drawnow;clear;close all
       
  % Compare test values to real values.
        
  function testLGCP
    values.real = load('realValuesLgcp.mat', 'x');
    values.test = load(strrep(which('test_lgcp.m'), 'test_lgcp.m', 'testValues/testLgcp.mat'), 'x');
    assertElementsAlmostEqual(mean(values.real.x), mean(values.test.x), 'relative', 0.05)

