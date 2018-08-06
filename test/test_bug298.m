function test_bug298

% MEM 1500mb
% WALLTIME 00:10:00

% TEST keyval keyvalx

warning('this test has become obsolete because ot the ft_getopt implementation');
return

keyVals  = {'frequency', 100,'order', 6,'decimation', 4,'steepness', 20,'name' ,'Stefan','matrix', magic(4)};
keyValsA = {'frequencyA',100,'orderA',6,'decimationA',4,'steepnessA',20,'nameA','Stefan','matrixA',magic(4)};
keyValsB = {'frequencyB',100,'orderB',6,'decimationB',4,'steepnessB',20,'nameB','Stefan','matrixB',magic(4)};
key = 'nameB';

kv = [keyVals keyValsA keyValsB];

disp('Timing of M-file and MEX-file without returning remaining key values:');
tic;
for k=1:1000;
  v = keyval(key, kv);
end
toc
tic;
for k=1:1000;
  vx = keyvalx(key, kv);
end
toc
        
disp('Timing of M-file and MEX-file *with* returning remaining key values:');        
        
tic;
for k=1:1000;
  [v,r] = keyval(key, kv);
end
toc
tic;
for k=1:1000;
  [vx,rx] = keyvalx(key, kv);
end
toc
                
