function tests = test_ft_connectivity_plm

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_plm

if nargout
  % assume that this is called by RUNTESTS
  tests = functiontests(localfunctions);
else
  % assume that this is called from the command line
  func = localfunctions;
  for i=1:numel(func)
    fprintf('evaluating %s\n', func2str(func{i}));
    feval(func{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_rpt_chan_time(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt   = 10;
nchan  = 3;
ntime  = 1000;

for i=1:nrpt
  input{i} = randn(nchan,ntime);
end

result = {};
result{end+1} = ft_connectivity_plm(input, 'fsample', 1000, 'bandwidth', 5);
result{end+1} = ft_connectivity_plm(input, 'fsample', 2000, 'bandwidth', 5);
result{end+1} = ft_connectivity_plm(input, 'fsample', 3000, 'bandwidth', 5);
result{end+1} = ft_connectivity_plm(input, 'fsample', 1000, 'bandwidth', 10);
result{end+1} = ft_connectivity_plm(input, 'fsample', 1000, 'bandwidth', 20);
result{end+1} = ft_connectivity_plm(input, 'fsample', 1000, 'bandwidth', 30);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

