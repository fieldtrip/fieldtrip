function tests = test_ft_specest_tfr

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_specest_tfr

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
verbose = false;
freqoi  = 1:500;

result = {};

result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 1, 'gwidth', 2);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 3, 'gwidth', 2);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 5, 'gwidth', 2);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 7, 'gwidth', 2);

result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 1, 'gwidth', 3);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 3, 'gwidth', 3);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 5, 'gwidth', 3);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 7, 'gwidth', 3);

result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 5, 'gwidth', 4, 'polyorder', 0);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 5, 'gwidth', 4, 'polyorder', 1);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 5, 'gwidth', 4, 'polyorder', 2);
result{end+1} = ft_specest_tfr(dat, time, 'verbose', verbose, 'timeoi', time, 'freqoi', freqoi, 'width', 5, 'gwidth', 4, 'polyorder', 3);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
