function inspect_qsubcellfun3

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_qsubcellfun3
% TEST qsubcellfun qsubfeval qsubget

% this should not run in the automated batch, because the torque queue
% will be completely full with other jobs, causing this job to timeout

if isempty(which('qsubcellfun'))
  [ftver, ftpath] = ft_version;
  addpath(fullfile(ftpath, 'qsub'));
end

result1 = cellfun(@subfunction, {1, 2, 3}, 'UniformOutput', false);

result2 = qsubcellfun(@subfunction, {1, 2, 3}, 'memreq', 1e8, 'timreq', 300, 'backend', 'local');
assert(isequal(result1, result2));

% % the following does not work, which is the correct behaviour
% % the subfunction cannot be located if passed as a string
% result2 = qsubcellfun('subfunction', {1, 2, 3}, 'memreq', 1e8, 'timreq', 300, 'backend', 'local');
% assert(isequal(result1, result2));

% this section was confirmed to work on 14 October 2012
result3 = qsubcellfun(@subfunction, {1, 2, 3}, 'memreq', 1e8, 'timreq', 300);
assert(isequal(result1, result3));

% this section fails on 14 October 2012
% result4 = qsubcellfun(@subsubfunction, {1, 2, 3}, 'memreq', 1e8, 'timreq', 300);
% assert(isequal(result1, result4));

function y = subfunction(x)
y = x.^2;

function y = subsubfunction(x)
y = subfunction(x);

