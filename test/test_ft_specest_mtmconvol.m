function tests = test_ft_specest_mtmconvol

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_specest_mtmconvol

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
freqoi  = 4:4:30;
tapsmofrq = 4*ones(size(freqoi));
timeoi  = time(5:10:end);
timwin  = 0.25 * ones(size(freqoi));
verbose = false;

result = {};

% vary the freqio and timwin
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi(2:2:end), 'timeoi', timeoi, 'timwin', timwin(2:2:end), 'taper', 'triang', 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi(3:3:end), 'timeoi', timeoi, 'timwin', timwin(3:3:end), 'taper', 'triang', 'verbose', verbose);

% vary the timeoi
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi(2:2:end), 'timwin', timwin, 'taper', 'triang', 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi(3:3:end), 'timwin', timwin, 'taper', 'triang', 'verbose', verbose);

% vary the taper
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'bartlett'      , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'barthannwin'   , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'blackman'      , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'blackmanharris', 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'bohmanwin'     , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'chebwin'       , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'flattopwin'    , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'gausswin'      , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hann'          , 'verbose', verbose); % hann and hanning are nearly the same in the MATLAB implementation, and exactly the same in external/signal
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'kaiser'        , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'nuttallwin'    , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'parzenwin'     , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'rectwin'       , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'taylorwin'     , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'tukeywin'      , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'triang'        , 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'alpha'         , 'verbose', verbose);

% vary the multitapering
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'tapsmofrq', 1*tapsmofrq, 'taper', 'dpss', 'verbose', verbose); % Warning: using only one taper for specified smoothing
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'tapsmofrq', 2*tapsmofrq, 'taper', 'dpss', 'verbose', verbose);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'tapsmofrq', 1*tapsmofrq, 'taper', 'sine', 'verbose', verbose); % Warning: using only one taper for specified smoothing
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'tapsmofrq', 2*tapsmofrq, 'taper', 'sine', 'verbose', verbose);

% vary the padding
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hanning', 'verbose', verbose, 'pad', 2, 'padtype', 'zero');
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hanning', 'verbose', verbose, 'pad', 3, 'padtype', 'zero');
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hanning', 'verbose', verbose, 'pad', 4, 'padtype', 'zero');

% try out the various padding types, don't remove the polynomial fit
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'zero');
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'mean');
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'localmean');
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'edge');
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hanning', 'verbose', verbose, 'polyorder', -1, 'pad', 1/0.9, 'padtype', 'mirror');

% the default polyorder 0 is already covered earlier
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hann', 'verbose', verbose, 'polyorder', 1);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hann', 'verbose', verbose, 'polyorder', 2);
result{end+1} = ft_specest_mtmconvol(dat, time, 'freqoi', freqoi, 'timeoi', timeoi, 'timwin', timwin, 'taper', 'hann', 'verbose', verbose, 'polyorder', 3);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
