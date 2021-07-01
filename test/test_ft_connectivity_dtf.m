function tests = test_ft_connectivity_dtf

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_dtf

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
function testOptions(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt  = 10;
nchan = 3;
nfreq = 4;
ntime = 5;

% construct something that mimics the transfer function
input = 10 * randn(nrpt, nchan, nchan, nfreq, ntime) + 1i * 10 * randn(nrpt, nchan, nchan, nfreq, ntime);
crsspctrm = zeros(nrpt, nchan, nchan, nfreq, ntime); % needed to test invfun

for rpt=1:nrpt
  for time=1:ntime
    for freq=1:nfreq
      fdat = 10 * randn(nchan, 1) + 1i * 10 * randn(nchan, 1);
      crsspctrm(rpt,:,:,freq,time) = fdat * ctranspose(fdat); % compute the cross-spectral density
    end
  end
end

result = {};
result{end+1} = ft_connectivity_dtf(input);
result{end+1} = ft_connectivity_dtf(input, 'crsspctrm', crsspctrm, 'invfun', 'inv'); % Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
result{end+1} = ft_connectivity_dtf(input, 'crsspctrm', crsspctrm, 'invfun', 'pinv');

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
