function test_qsubcellfun3

% TEST test_qsubcellfun3
% TEST qsubcellfun qsubfeval qsubget


result1 = cellfun(@subfunction, {1, 2, 3}, 'UniformOutput', false);

result2 = qsubcellfun(@subfunction, {1, 2, 3}, 'memreq', 1e8, 'timreq', 300, 'backend', 'local');
assert(isequal(result1, result2));

% the following does not work, which is the correct behaviour
% the subfunction cannot be located if passed as a string
% result2 = qsubcellfun('subfunction', {1, 2, 3}, 'memreq', 1e8, 'timreq', 300, 'backend', 'local');
% assert(isequal(result1, result2));

if false
  % this should not run in the automated batch, because the torque queue
  % will be completely full with other jobs, causing this job to timeout
  
  % this section was confirmed to work on 14 October 2012
  result3 = qsubcellfun(@subfunction, {1, 2, 3}, 'memreq', 1e8, 'timreq', 300);
  assert(isequal(result1, result3));
end

if false
  % this should not run in the automated batch, because the torque queue
  % will be completely full with other jobs, causing this job to timeout
  
  % this section fails on 14 October 2012
  result4 = qsubcellfun(@subsubfunction, {1, 2, 3}, 'memreq', 1e8, 'timreq', 300);
  assert(isequal(result1, result4));
end


function y = subfunction(x)
y = x.^2;

function y = subsubfunction(x)
y = subfun(x);

