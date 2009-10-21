function [freq] = freqanalysis_mtmwelch(cfg, data);

% FREQANALYSIS_MTMWELCH performs frequency analysis on any time series
% trial data using the 'multitaper method' (MTM) based on discrete
% prolate spheroidal sequences (Slepian sequences) as tapers. Alternatively,
% you can use conventional tapers (e.g. Hanning).
% 
% Besides multitapering, this function uses Welch's averaged, modified
% periodogram method. The data is divided into a number of sections with
% overlap, each section is windowed with the specified taper(s) and the
% powerspectra are computed and averaged over the sections in each trial.
%
% Use as
%   [freq] = freqanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from
% the PREPROCESSING function. The configuration should be according to
%   cfg.method     = method used for frequency or time-frequency decomposition
%                    see FREQANALYSIS for details
%   cfg.output     = 'pow'       return the power-spectra
%                    'powandcsd' return the power and the cross-spectra
%   cfg.taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%
% For cfg.output='powandcsd', you should specify the channel combinations
% between which to compute the cross-spectra as cfg.channelcmb. Otherwise
% you should specify only the channels in cfg.channel.
% 
%   cfg.channel    = Nx1 cell-array with selection of channels (default = 'all'),
%                    see CHANNELSELECTION for details
%   cfg.channelcmb = Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                    see CHANNELCOMBINATION for details
%
% This function uses FREQANALYSIS_MTMCONVOL for the low-level
% computations, and you can use the options of that function to specify
% the length of the time windows, the amount of overlap, and the amount
% of spectral smoothing (in case of dpss tapers) per window.
%
%   cfg.foi        = vector 1 x numfoi, frequencies of interest
%   cfg.t_ftimwin  = vector 1 x numfoi, length of time window (in seconds)
%   cfg.tapsmofrq  = vector 1 x numfoi, the amount of spectral smoothing through
%                    multi-tapering. Note that 4 Hz smoothing means
%                    plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
%   cfg.toi        = vector 1 x numtoi, the times on which the analysis windows 
%                    should be centered (in seconds)
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrials = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.keeptapers = 'yes' or 'no', return individual tapers or average (default = 'no')
%   cfg.pad        = number or 'maxperlen', length in seconds to which the data can be padded out (default = 'maxperlen')
%
% See also FREQANALYSIS_MTMCONVOL, FREQANALYSIS

% This function depends on FREQANALYSIS which uses cfg.method = 'mtmconvol'

% Copyright (C) 2005-2006, F.C. Donders Centre
%
% $Log: freqanalysis_mtmwelch.m,v $
% Revision 1.12  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.11  2008/01/18 13:14:50  sashae
% added option for trial selection, updated documentation
%
% Revision 1.10  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.9  2006/06/20 16:25:58  ingnie
% updated documentation
%
% Revision 1.8  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.7  2006/05/23 16:05:20  ingnie
% updated documentation
%
% Revision 1.6  2006/05/04 14:26:35  roboos
% added some real documentation
%
% Revision 1.5  2006/03/14 08:09:22  roboos
% added copyrigth and cvs log statement
% 

fieldtripdefs

% ensure that this function is started as a subfunction of the FREQANALYSIS wrapper
if ~exist('OCTAVE_VERSION')
  [s, i] = dbstack;
  if length(s)>1
    [caller_path, caller_name, caller_ext] = fileparts(s(2).name);
  else
    caller_path = '';
    caller_name = '';
    caller_ext  = '';
  end
  % evalin('caller', 'mfilename') does not work for Matlab 6.1 and 6.5
  if ~strcmp(caller_name, 'freqanalysis')
    error(['you should call FREQANALYSIS, instead of the ' upper(mfilename) ' subfunction']);
  end
end

tmin =  inf;
tmax = -inf;
for i=1:length(data.time)
  tmin = min(tmin, data.time{i}(1));
  tmax = max(tmax, data.time{i}(end));
end

fprintf('taking every sample as time of interest\n');

cfgconvol        = cfg;
cfgconvol.method = 'mtmconvol';
cfgconvol.toi    = tmin:(1/data.fsample):tmax;
cfgconvol.trials = 'all'; % trial selection already applied during first call of freqanalysis

% use mtmconvol to do the dirty work
freq = freqanalysis(cfgconvol, data);

% determine the time dimension
if strcmp(freq.dimord, 'chan_freq_time')
  timedim = 3;
elseif strcmp(freq.dimord, 'rpt_chan_freq_time')
  timedim = 4;
else
  error('unexpected dimord');
end

% average over time
% NOTE: the degrees of freedom should be returned by freqanalysis_mtmconvol, 
% and used here to weigh every trial and timepoint accordingly. But
% currently that is not yet possible.
freq.powspctrm = nan_mean(freq.powspctrm, timedim);
if isfield(freq, 'crsspctrm')
  freq.crsspctrm = nan_mean(freq.crsspctrm, timedim);
end

% remove the time axis
freq = rmfield(freq, 'time');

% update the dimord
if strcmp(freq.dimord, 'chan_freq_time')
  freq.dimord = 'chan_freq';
elseif strcmp(freq.dimord, 'rpt_chan_freq_time')
  freq.dimord = 'rpt_chan_freq';
end

% only update the method and trials fields
freq.cfg.method = cfg.method;
freq.cfg.trials = cfg.trials;
