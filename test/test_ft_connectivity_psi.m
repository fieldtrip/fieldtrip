function tests = test_ft_connectivity_psi

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_psi
% DATA no

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

nrpt = 10;
nchan = 2;
nfreq = 501;
fftdata = zeros(nrpt, nchan, nfreq);
csddata = zeros(nrpt, nchan, nchan, nfreq);
for i = 1:nrpt
  % simulate data as two time shifted copies
  x = randn(1,1010);
  y = randn(2,1000);
  y(1,:) = x(1:1000)  + y(1,:)./10;
  y(2,:) = x(11:1010) + y(2,:)./10;
  f = fft(y, [], 2);
  fftdata(i,:,:) = f(:, 1:501);
  csddata(i,1,1,:) = f(1, 1:501).*conj(f(1, 1:501));
  csddata(i,1,2,:) = f(1, 1:501).*conj(f(2, 1:501));
  csddata(i,2,1,:) = f(2, 1:501).*conj(f(1, 1:501));
  csddata(i,2,2,:) = f(2, 1:501).*conj(f(2, 1:501));
end
dimord  = 'rpt_chan_chan_freq';

result = {};
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'nbin', []);
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'nbin', 2);
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'nbin', 4);
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'nbin', [], 'normalize', 'yes');
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'nbin', 2,  'normalize', 'yes');
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'nbin', 4,  'normalize', 'yes');

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    if all(result{i}(:)==0) && all(result{j}(:)==0)
      continue;
    end
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_rpt_chancmb(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrpt = 10;
nchan = 2;
nfreq = 501;
fftdata = zeros(nrpt, nchan,   nfreq);
csddata = zeros(nrpt, nchan+1, nfreq);
for i = 1:nrpt
  x = randn(1,1010);
  y = randn(2,1000);
  y(1,:) = x(1:1000)  + y(1,:)./10;
  y(2,:) = x(11:1010) + y(2,:)./10;
  f = fft(y, [], 2);
  fftdata(i,:,:) = f(:, 1:501);
end
csddata(:,1,:) = fftdata(:,1,:).*conj(fftdata(:,2,:));
csddata(:,2,:) = fftdata(:,1,:).*conj(fftdata(:,1,:));
csddata(:,3,:) = fftdata(:,2,:).*conj(fftdata(:,2,:));

powindx = [2 3;2 2;3 3];
dimord  = 'rpt_chancmb_freq';

result = {};
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'powindx', powindx, 'nbin', []);
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'powindx', powindx, 'nbin', 2);
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'powindx', powindx, 'nbin', 4);
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'powindx', powindx, 'nbin', [], 'normalize', 'yes');
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'powindx', powindx, 'nbin', 2,  'normalize', 'yes');
result{end+1} = ft_connectivity_psi(csddata, 'dimord', dimord, 'powindx', powindx, 'nbin', 4,  'normalize', 'yes');

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    if all(result{i}(:)==0) && all(result{j}(:)==0)
      continue;
    end
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

