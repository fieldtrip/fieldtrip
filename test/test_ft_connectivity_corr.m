function tests = test_ft_connectivity_corr

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_corr

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

nrpt   = 10;
nchan  = 3;
ntime  = 100;
dimord = 'rpt_chan_chan';

input = zeros(nrpt, nchan, nchan);
for i=1:nrpt
  dat = 10 * randn(nchan, ntime);
  input(i,:,:) = cov(dat');
end

result = {};
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', false);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 1);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 1:2);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 1:3);

assert(all(diag(result{2})==1));
assert(all(size(result{2})==[3 3]));
assert(all(size(result{3})==[2 2]));
assert(all(size(result{4})==[1 1]));
assert(all(size(result{5})==[0 0]));

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_rpt_chan_chan_time(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt   = 10;
nchan  = 3;
nseg   = 20;  % sliding window, there are 20 segments per trial
ntime  = 50; % sliding window, there are 50 samples per segment
dimord = 'rpt_chan_chan_time';

input = zeros(nrpt, nchan, nchan);
for rpt=1:nrpt
  for seg=1:nseg
    dat = 10 * randn(nchan, ntime);
    input(rpt,:,:,seg) = cov(dat');
  end
end


result = {};
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', false);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 1, 'allchanindx', 1:3);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 2, 'allchanindx', 1:3);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 1:2, 'allchanindx', 1:3);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_rpt_chan_chan_freq(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt   = 10;
nchan  = 3;
nfreq  = 5;
dimord = 'rpt_chan_chan_freq';


input = zeros(nrpt, nchan, nchan, nfreq);
for rpt=1:nrpt
  for freq=1:nfreq
    fdat = 10 * randn(nchan, 1) + 1i * 10 * randn(nchan, 1);
    input(rpt,:,:,freq) = fdat * ctranspose(fdat); % compute the cross-spectral density
  end
end

result = {};
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', false);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 1, 'allchanindx', 1:3);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 2, 'allchanindx', 1:3);
result{end+1} = ft_connectivity_corr(input, 'dimord', dimord, 'hasjack', false, 'pownorm', true, 'pchanindx', 1:2, 'allchanindx', 1:3);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
