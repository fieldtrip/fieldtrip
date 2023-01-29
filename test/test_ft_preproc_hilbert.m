function tests = test_ft_preproc_hilbert

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_hilbert

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

handlenan = [];
padnan = [];

result = [];
result{end+1} = ft_preproc_hilbert(dat, 'abs'     , handlenan, padnan);
result{end+1} = ft_preproc_hilbert(dat, 'complex' , handlenan, padnan);
result{end+1} = ft_preproc_hilbert(dat, 'real'    , handlenan, padnan);
result{end+1} = ft_preproc_hilbert(dat, 'imag'    , handlenan, padnan);
result{end+1} = ft_preproc_hilbert(dat, 'absreal' , handlenan, padnan);
result{end+1} = ft_preproc_hilbert(dat, 'absimag' , handlenan, padnan);
result{end+1} = ft_preproc_hilbert(dat, 'angle'   , handlenan, padnan);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end