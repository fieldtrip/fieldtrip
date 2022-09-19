function tests = test_ft_preproc_rereference

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_rereference

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

handlenan = false;
leadfield = randn(8, 1000);

result = {};
result{end+1} = ft_preproc_rereference(dat, 1,   'avg', handlenan);
result{end+1} = ft_preproc_rereference(dat, 1:4, 'avg', handlenan);
result{end+1} = ft_preproc_rereference(dat, 1:8, 'avg', handlenan);
result{end+1} = ft_preproc_rereference(dat, 1:4, 'median', handlenan);
result{end+1} = ft_preproc_rereference(dat, 1:8, 'median', handlenan);
result{end+1} = ft_preproc_rereference(dat, 1,   'rest', handlenan, leadfield);
result{end+1} = ft_preproc_rereference(dat, 1:4, 'rest', handlenan, leadfield);
result{end+1} = ft_preproc_rereference(dat, 1:8, 'rest', handlenan, leadfield);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end