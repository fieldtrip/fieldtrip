function tests = test_ft_preproc_online_downsample

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_online_downsample_init ft_preproc_online_downsample_apply

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

state1 = ft_preproc_online_downsample_init(1); % factor 1x
state2 = ft_preproc_online_downsample_init(2); % factor 2x
state3 = ft_preproc_online_downsample_init(3); % factor 3x

result = {};
[state1, result{end+1}] = ft_preproc_online_downsample_apply(state1, dat); % this will be the same on each call

[state2, result{end+1}] = ft_preproc_online_downsample_apply(state2, dat); % this will be the same on each call

[state3, result{end+1}] = ft_preproc_online_downsample_apply(state3, dat);
[state3, result{end+1}] = ft_preproc_online_downsample_apply(state3, dat);
[state3, result{end+1}] = ft_preproc_online_downsample_apply(state3, dat);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
