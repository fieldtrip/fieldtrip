function tests = test_ft_preproc_bandstopfilter

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_bandstopfilter

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

dat = randn(nchan, nsample) + 1;
Fs  = 1000;
Fbs = [48 52];

instabilityfix  = [];
df              = []; 
wintype         = []; 
dev             = []; 
plotfiltresp    = [];
usefftfilt      = [];

result = {};
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, [], 'brickwall', [], instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt); % this does not use order and dir
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 1, 'but'      , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 1, 'but'      , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 1, 'but'      , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 1, 'but'      , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 1, 'but'      , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'but'      , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'but'      , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'but'      , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'but'      , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'but'      , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firws'    , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firws'    , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firws'    , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firws'    , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firws'    , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firws'    , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firws'    , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firws'    , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'fir'      , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'fir'      , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'fir'      , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'fir'      , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'fir'      , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'fir'      , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'fir'      , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'fir'      , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firls'    , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firls'    , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firls'    , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firls'    , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firls'    , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firls'    , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firls'    , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 2, 'firls'    , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'but'      , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'but'      , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'but'      , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'but'      , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'but'      , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firws'    , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firws'    , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firws'    , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firws'    , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firws'    , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firws'    , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firws'    , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firws'    , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'fir'      , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'fir'      , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'fir'      , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'fir'      , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'fir'      , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'fir'      , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'fir'      , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'fir'      , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firls'    , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firls'    , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firls'    , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firls'    , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firls'    , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firls'    , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firls'    , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 4, 'firls'    , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'but'      , 'onepass'                  , 'split', [], [], [], [], []); % this one is instable
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'but'      , 'onepass-reverse'          , 'split', [], [], [], [], []);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'but'      , 'twopass'                  , 'split', [], [], [], [], []);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'but'      , 'twopass-reverse'          , 'split', [], [], [], [], []);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'but'      , 'twopass-average'          , 'split', [], [], [], [], []);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firws'    , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firws'    , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firws'    , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firws'    , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firws'    , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firws'    , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firws'    , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firws'    , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'fir'      , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'fir'      , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'fir'      , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'fir'      , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'fir'      , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'fir'      , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'fir'      , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'fir'      , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firls'    , 'onepass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firls'    , 'onepass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firls'    , 'onepass-zerophase'        , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firls'    , 'onepass-reverse-zerophase', instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firls'    , 'onepass-minphase'         , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firls'    , 'twopass'                  , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firls'    , 'twopass-reverse'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);
result{end+1} = ft_preproc_bandstopfilter(dat, Fs, Fbs, 8, 'firls'    , 'twopass-average'          , instabilityfix, df, wintype, dev, plotfiltresp, usefftfilt);

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequal(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end
