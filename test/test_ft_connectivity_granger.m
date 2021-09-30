function tests = test_ft_connectivity_granger

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_granger

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

nrpt = 10;
nchan = 4;
nfreq = 5;

H = randn(nrpt, nchan, nchan, nfreq) + 1i * randn(nrpt, nchan, nchan, nfreq);
Z = randn(nrpt, nchan, nchan);
S = randn(nrpt, nchan, nchan, nfreq) + 1i * randn(nrpt, nchan, nchan, nfreq);
dimord = 'rpt_chan_chan_freq';

result = {};
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'granger'      , 'hasjack', false);
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'instantaneous', 'hasjack', false);
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'total'        , 'hasjack', false);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_rpt_chan_chan_freq_time(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt = 10;
nchan = 4;
nfreq = 5;
ntime = 6;

H = randn(nrpt, nchan, nchan, nfreq, ntime) + 1i * randn(nrpt, nchan, nchan, nfreq, ntime);
Z = randn(nrpt, nchan, nchan, ntime);
S = randn(nrpt, nchan, nchan, nfreq, ntime) + 1i * randn(nrpt, nchan, nchan, nfreq, ntime);
dimord = 'rpt_chan_chan_freq_time';

result = {};
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'granger'      , 'hasjack', false);
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'instantaneous', 'hasjack', false);
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'total'        , 'hasjack', false);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_rpt_chan_freq(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt = 10;
nchan = 4;
nfreq = 5;

H = randn(nrpt, nchan*nchan, nfreq) + 1i * randn(nrpt, nchan*nchan, nfreq);
Z = randn(nrpt, nchan*nchan);
S = randn(nrpt, nchan*nchan, nfreq) + 1i * randn(nrpt, nchan*nchan, nfreq);
dimord = 'rpt_chan_freq_time';

result = {};
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'granger'      , 'hasjack', false);
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'instantaneous', 'hasjack', false);
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'total'        , 'hasjack', false);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_rpt_chan_freq_time(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt = 10;
nchan = 4;
nfreq = 5;
ntime = 6;

H = randn(nrpt, nchan*nchan, nfreq, ntime) + 1i * randn(nrpt, nchan*nchan, nfreq, ntime);
Z = randn(nrpt, nchan*nchan, ntime);
S = randn(nrpt, nchan*nchan, nfreq, ntime) + 1i * randn(nrpt, nchan*nchan, nfreq, ntime);
dimord = 'rpt_chan_freq_time';

result = {};
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'granger'      , 'hasjack', false);
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'instantaneous', 'hasjack', false);
result{end+1} = ft_connectivity_granger(H, Z, S, 'dimord', dimord, 'method', 'total'        , 'hasjack', false);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
