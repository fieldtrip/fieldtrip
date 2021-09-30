function tests = test_ft_preproc_online_filter

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_online_filter_init ft_preproc_online_filter_apply

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

nchan   = 2;
nsample = 100;

result = {};

[B,A] = butter(6, 0.05, 'high');
state = ft_preproc_online_filter_init(B, A, zeros(nchan,1));

[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);

[B,A] = butter(6, 0.20, 'low');
state = ft_preproc_online_filter_init(B, A, zeros(nchan,1));

[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);

[B,A] = butter(6, [0.05 0.20]);
state = ft_preproc_online_filter_init(B, A, zeros(nchan,1));

[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);
[state, result{end+1}] = ft_preproc_online_filter_apply(state, randn(nchan, nsample) + 10);

% the first  part is high-pass filtered
% the second part is low-pass  filtered
% the third  part is band-pass filtered
plot(cat(2, result{:})');

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
