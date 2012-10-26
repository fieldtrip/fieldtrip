% Run all test scripts in this directory and report results

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: run_tests.m,v 1.7 2005/08/17 09:03:12 kosti Exp $

files = dir('.');
count_total = 0;
count_failed = 0;
for i = 1:length(files)
    file = files(i).name;
    % for ML6.5+: if regexp(file,'^test_.*\.m$')
    l = length(file);
    if l>8
      if strcmp(file([1:5,l-1:l]), 'test_.m')
          count_total = count_total + 1;
          fprintf('Test %s',file);
          file = file(1:length(file)-2);
          try
              T=evalc(file);
              %fprintf('%s',T);
              fprintf(' PASSED\n', file);
          catch
              fprintf(' FAILED, reason:\n  %s\n', lasterr);
              count_failed=count_failed+1;
          end
      end
    end
end

fprintf('----\nTest results:\n');
fprintf('  Number of tests %d\n', count_total);
fprintf('  PASSED          %d\n', count_total-count_failed);
fprintf('  FAILED          %d\n', count_failed);

