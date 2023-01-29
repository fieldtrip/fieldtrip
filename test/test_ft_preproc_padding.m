function tests = test_ft_preproc_padding

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_padding

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

result = [];
result{end+1} = ft_preproc_padding(dat, 'zero',      10);
result{end+1} = ft_preproc_padding(dat, 'mean',      10);
result{end+1} = ft_preproc_padding(dat, 'localmean', 10);
result{end+1} = ft_preproc_padding(dat, 'edge',      10);
result{end+1} = ft_preproc_padding(dat, 'mirror',    10);
result{end+1} = ft_preproc_padding(dat, 'nan',       10);
result{end+1} = ft_preproc_padding(dat, 'remove',    10);
result{end+1} = ft_preproc_padding(dat, 'zero',      20);
result{end+1} = ft_preproc_padding(dat, 'mean',      20);
result{end+1} = ft_preproc_padding(dat, 'localmean', 20);
result{end+1} = ft_preproc_padding(dat, 'edge',      20);
result{end+1} = ft_preproc_padding(dat, 'mirror',    20);
result{end+1} = ft_preproc_padding(dat, 'nan',       20);
result{end+1} = ft_preproc_padding(dat, 'remove',    20);
result{end+1} = ft_preproc_padding(dat, 'zero',      10, 20);
result{end+1} = ft_preproc_padding(dat, 'mean',      10, 20);
result{end+1} = ft_preproc_padding(dat, 'localmean', 10, 20);
result{end+1} = ft_preproc_padding(dat, 'edge',      10, 20);
result{end+1} = ft_preproc_padding(dat, 'mirror',    10, 20);
result{end+1} = ft_preproc_padding(dat, 'nan',       10, 20);
result{end+1} = ft_preproc_padding(dat, 'remove',    10, 20);

for i=1:numel(result)
  % the number of samples should be larger, or smaller upon 'remove'
  assert(size(result{i},2)~=nsample);
end

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
