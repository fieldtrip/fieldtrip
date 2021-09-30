function tests = test_ft_connectivity_pdc

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_pdc

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
function test_rpt_chan_chan_freq(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt   = 10;
nchan  = 3;
nfreq  = 4;
H      = randn(nrpt, nchan, nchan, nfreq);

dat = randn(nchan,1000);
noisecov = cov(dat');

result = {};
result{end+1} = ft_connectivity_pdc(H, 'noisecov', []);
result{end+1} = ft_connectivity_pdc(H, 'noisecov', noisecov);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_rpt_chan_chan_freq_time(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt   = 10;
nchan  = 3;
nfreq  = 4;
ntime  = 5;
H      = randn(nrpt, nchan, nchan, nfreq, ntime);

dat = randn(nchan,1000);
noisecov = cov(dat');

result = {};
result{end+1} = ft_connectivity_pdc(H, 'noisecov', []);
result{end+1} = ft_connectivity_pdc(H, 'noisecov', noisecov);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
