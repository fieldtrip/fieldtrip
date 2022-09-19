function tests = test_ft_specest_hilbert

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_specest_hilbert

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
freqoi  = (10:10:100)';
width   = 5;

result = {};

% first explore the different filter types and orders
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but',   'bpfiltord', 2,   'bpfiltdir', 'twopass', 'timeoi', time);
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but',   'bpfiltord', 2,   'bpfiltdir', 'twopass', 'timeoi', time(1:2:end));
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but',   'bpfiltord', 3,   'bpfiltdir', 'twopass', 'timeoi', time(1:2:end));

% bpfiltdir=twopass is the default, use a different filter order for the following ones
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 1, 'bpfiltdir', 'onepass');
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 1, 'bpfiltdir', 'onepass-reverse');
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 1, 'bpfiltdir', 'twopass');
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 1, 'bpfiltdir', 'twopass-reverse');
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 1, 'bpfiltdir', 'twopass-average');

% padtype=zero is the default, use a different filter order for the following ones
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 2, 'bpfiltdir', 'twopass', 'pad', 2, 'padtype', 'zero'      );
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 2, 'bpfiltdir', 'twopass', 'pad', 2, 'padtype', 'mean'      );
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 2, 'bpfiltdir', 'twopass', 'pad', 2, 'padtype', 'localmean' );
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 2, 'bpfiltdir', 'twopass', 'pad', 2, 'padtype', 'edge'      );
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 2, 'bpfiltdir', 'twopass', 'pad', 2, 'padtype', 'mirror'    );

% polyorder=0 is the default, use a different filter order for the following ones
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 3, 'bpfiltdir', 'twopass', 'polyorder', 0);
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 3, 'bpfiltdir', 'twopass', 'polyorder', 1);
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 3, 'bpfiltdir', 'twopass', 'polyorder', 2);
result{end+1} = ft_specest_hilbert(dat, time, 'verbose', verbose, 'freqoi', freqoi, 'width', width, 'bpfilttype', 'but', 'bpfiltord', 3, 'bpfiltdir', 'twopass', 'polyorder', 3);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
