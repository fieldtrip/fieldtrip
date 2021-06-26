function tests = test_ft_preproc_resample

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_resample

if nargout
  % assume that this is called by RUNTESTS
  tests = functiontests(localfunctions);
else
  % assume that this is called from the command line
  fn = localfunctions;
  for i=1:numel(fn)
    feval(fn{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testOptions(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nchan   = 8;
nsample = 1000;
dat     = randn(nchan, nsample) + 1;

result = {};
result{end+1} = ft_preproc_resample(dat, 1000, 100, 'resample');
result{end+1} = ft_preproc_resample(dat, 1000, 250, 'resample');
result{end+1} = ft_preproc_resample(dat, 1000, 500, 'resample');
result{end+1} = ft_preproc_resample(dat, 1000, 100, 'decimate');
result{end+1} = ft_preproc_resample(dat, 1000, 250, 'decimate');
result{end+1} = ft_preproc_resample(dat, 1000, 500, 'decimate');
result{end+1} = ft_preproc_resample(dat, 1000, 100, 'downsample');
result{end+1} = ft_preproc_resample(dat, 1000, 250, 'downsample');
result{end+1} = ft_preproc_resample(dat, 1000, 500, 'downsample');
result{end+1} = ft_preproc_resample(dat, 1000, 100, 'fft');
result{end+1} = ft_preproc_resample(dat, 1000, 250, 'fft');
result{end+1} = ft_preproc_resample(dat, 1000, 500, 'fft');

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end