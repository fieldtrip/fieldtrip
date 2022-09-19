function tests = test_ft_preproc_denoise

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_denoise

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

nsample = 1000;
dat     = randn(8, nsample) + 1;
refdat  = randn(4, nsample) + 1;

result = [];
result{end+1} = ft_preproc_denoise(dat, refdat(1,:), false);
result{end+1} = ft_preproc_denoise(dat, refdat(2,:), false);
result{end+1} = ft_preproc_denoise(dat, refdat(3,:), false);
result{end+1} = ft_preproc_denoise(dat, refdat(4,:), false);
result{end+1} = ft_preproc_denoise(dat, refdat(1,:), true);
result{end+1} = ft_preproc_denoise(dat, refdat(2,:), true);
result{end+1} = ft_preproc_denoise(dat, refdat(3,:), true);
result{end+1} = ft_preproc_denoise(dat, refdat(4,:), true);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end