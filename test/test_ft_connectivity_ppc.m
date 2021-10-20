function tests = test_ft_connectivity_ppc

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_ppc

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
function test_rpt_chan_chan(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the code does exactly the same thing for all supported dimords, so no need to test all of them
%   rpt_chan_chan
%   rpt_chan_chan_freq
%   rpt_chan_chan_freq_time
%   rpt_chancmb
%   rpt_chancmb_freq
%   rpt_chancmb_freq_time

nrpt  = 10;
nchan = 3;
input = randn(nrpt, nchan, nchan) + 1i*randn(nrpt, nchan, nchan);

result = {};
result{end+1} = ft_connectivity_ppc(input, 'weighted', true);
result{end+1} = ft_connectivity_ppc(input, 'weighted', false);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
