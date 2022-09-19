function tests = test_ft_specest_irasa

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_specest_irasa

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
nchan   = 2;
nsample = 1000;

dat     = randn(nchan, nsample) + 1 + linspace(0,1,nsample);
time    = (0:(nsample-1))/fsample;
freqoi  = 1:150;

result = {};
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'original');
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal');
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi(2:2:end), 'output', 'original');
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi(2:2:end), 'output', 'fractal');

% try out another padding
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'original', 'pad', 2);
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal',  'pad', 2);

% try out the various padding types, don't remove the polynomial fit
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal', 'polyorder', -1, 'pad', 2, 'padtype', 'zero');
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal', 'polyorder', -1, 'pad', 2, 'padtype', 'mean');
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal', 'polyorder', -1, 'pad', 2, 'padtype', 'localmean');
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal', 'polyorder', -1, 'pad', 2, 'padtype', 'edge');
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal', 'polyorder', -1, 'pad', 2, 'padtype', 'mirror');

% the default polyorder 1 is already covered earlier
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal', 'polyorder', 0);
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal', 'polyorder', 2);
result{end+1} = ft_specest_irasa(dat, time, 'verbose', false, 'freqoi', freqoi, 'output', 'fractal', 'polyorder', 3);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
