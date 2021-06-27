function tests = test_ft_specest_neuvar

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_specest_neuvar

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

fsample = 1000;
nchan   = 8;
nsample = 1000;

dat     = randn(nchan, nsample) + 1 + linspace(0,1,nsample);
time    = (0:(nsample-1))/fsample;

result = {};

result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false);
result{end+1} = ft_specest_neuvar(dat, time, 'order', 2, 'verbose', false);
result{end+1} = ft_specest_neuvar(dat, time, 'order', 3, 'verbose', false);

% the default polyorder 1 is already covered earlier
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'polyorder', 0);
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'polyorder', 2);
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'polyorder', 3);

% try out various amounts of padding
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'pad', 2, 'padtype', 'zero');
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'pad', 3, 'padtype', 'zero');
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'pad', 4, 'padtype', 'zero');

% try out the various padding types, don't remove the polynomial fit
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'polyorder', -1, 'pad', 2.5, 'padtype', 'zero');
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'polyorder', -1, 'pad', 2.5, 'padtype', 'mean');
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'polyorder', -1, 'pad', 2.5, 'padtype', 'localmean');
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'polyorder', -1, 'pad', 2.5, 'padtype', 'edge');
result{end+1} = ft_specest_neuvar(dat, time, 'order', 1, 'verbose', false, 'polyorder', -1, 'pad', 2.5, 'padtype', 'mirror');

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
