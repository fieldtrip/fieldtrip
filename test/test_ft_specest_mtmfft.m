function tests = test_ft_specest_mtmfft

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_specest_mtmfft

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
freqoi  = 1:500;
verbose = false;

result = {};
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'bartlett'      , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'barthannwin'   , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'blackman'      , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'blackmanharris', 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'bohmanwin'     , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'chebwin'       , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'flattopwin'    , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'gausswin'      , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hann'          , 'verbose', verbose); % hann and hanning are nearly the same in the MATLAB implementation, and exactly the same in external/signal
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'kaiser'        , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'nuttallwin'    , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'parzenwin'     , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'rectwin'       , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'taylorwin'     , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'tukeywin'      , 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'triang'        , 'verbose', verbose);

result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', 11:20, 'taper', 'hanning', 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', 21:30, 'taper', 'hanning', 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', 31:40, 'taper', 'hanning', 'verbose', verbose);

result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'tapsmofrq', 1, 'taper', 'dpss', 'verbose', verbose); % Warning: using only one taper for specified smoothing
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'tapsmofrq', 3, 'taper', 'dpss', 'verbose', verbose);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'tapsmofrq', 5, 'taper', 'dpss', 'verbose', verbose);

result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hanning', 'verbose', verbose, 'pad', 1.1, 'padtype', 'zero');
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hanning', 'verbose', verbose, 'pad', 1.2, 'padtype', 'zero');
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hanning', 'verbose', verbose, 'pad', 1.3, 'padtype', 'zero');

% try out the various padding types, don't remove the polynomial fit
% this will result in warning: output frequencies are different from input frequencies
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'zero');
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'mean');
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'localmean');
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'edge');
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', freqoi, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'mirror');

% the default polyorder 0 was already computed above, so use different frequencies
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', 2:2:100, 'taper', 'hanning', 'verbose', verbose, 'polyorder', 0);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', 2:2:100, 'taper', 'hanning', 'verbose', verbose, 'polyorder', 1);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', 2:2:100, 'taper', 'hanning', 'verbose', verbose, 'polyorder', 2);
result{end+1} = ft_specest_mtmfft(dat, time, 'freqoi', 2:2:100, 'taper', 'hanning', 'verbose', verbose, 'polyorder', 3);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
