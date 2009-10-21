function [freq] = freqanalysis(cfg, data);

% FREQANALYSIS performs frequency and time-frequency analysis
% on time series data over multiple trials
%
% Use as
%   [freq] = freqanalysis(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the PREPROCESSING function. The configuration depends on the type
% of computation that you want to perform.
%
% The configuration should contain:
%   cfg.method     = different methods of calculating the spectra
%                    'mtmfft', analyses an entire spectrum for the entire data
%                     length, implements multitaper frequency transformation
%                    'mtmconvol', implements multitaper time-frequency transformation
%                     based on multiplication in the frequency domain
%                    'mtmwelch', performs frequency analysis using Welch's averaged
%                     modified periodogram method of spectral estimation
%                    'wltconvol', implements wavelet time frequency transformation
%                     (using Morlet wavelets) based on multiplication in the frequency domain
%                    'tfr', implements wavelet time frequency transformation
%                     (using Morlet wavelets) based on convolution in the time domain
%
% The other cfg options depend on the method that you select. You should
% read the help of the respective subfunction FREQANALYSIS_XXX for the
% corresponding parameter options and for a detailed explanation of each method.
%
% See also FREQANALYSIS_MTMFFT, FREQANALYSIS_MTMCONVOL, FREQANALYSIS_MTMWELCH
% FREQANALYSIS_WLTCONVOL, FREQANALYSIS_TFR

% Undocumented local options:
% cfg.label
% cfg.labelcmb
% cfg.sgn
% cfg.sgncmb

% Copyright (C) 2003-2006, F.C. Donders Centre, Pascal Fries
% Copyright (C) 2004-2006, F.C. Donders Centre, Markus Siegel
%
% $Log: freqanalysis.m,v $
% Revision 1.46  2009/09/30 12:50:16  jansch
% added mvar as a supported datatype
%
% Revision 1.45  2009/07/30 19:44:43  ingnie
% also allow datatype comp (by chekdata)
%
% Revision 1.44  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.43  2009/01/12 13:05:20  sashae
% small change in call to checkconfig
%
% Revision 1.42  2008/09/23 12:30:33  sashae
% checkconfig: checks if the input cfg is valid for this function
%
% Revision 1.41  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.40  2008/05/06 15:43:46  sashae
% change in trial selection, cfg.trials can be logical
%
% Revision 1.39  2008/01/18 13:14:50  sashae
% added option for trial selection, updated documentation
%
% Revision 1.38  2007/11/05 09:48:56  roboos
% removed comments that pertain to the matlab compiler
%
% Revision 1.37  2007/05/29 16:58:06  ingnie
% moved making of offset field to checkdata by adding 'hasoffset' is 'yes'
%
% Revision 1.36  2007/05/02 15:59:13  roboos
% be more strict on the input and output data: It is now the task of
% the private/checkdata function to convert the input data to raw
% data (i.e. as if it were coming straight from preprocessing).
% Furthermore, the output data is NOT converted back any more to the
% input data, i.e. the output data is the same as what it would be
% on raw data as input, regardless of the actual input.
%
% Revision 1.35  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.34  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.33  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.32  2006/06/20 16:22:22  ingnie
% removed handling cfg.channel and cfg.channelcmb to freqanalysis_xxx functions
%
% Revision 1.31  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.30  2006/05/23 16:05:20  ingnie
% updated documentation
%
% Revision 1.29  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.28  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.27  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.26  2005/08/01 12:41:22  roboos
% added some pragma-like comments that help the matlab compiler with the dependencies
%
% Revision 1.25  2005/07/03 20:15:30  roboos
% replaced eval() by feval() for better support by Matlab compiler
%
% Revision 1.24  2005/06/02 12:24:47  roboos
% replaced the local conversion of average into raw trials by a call to the new helper function DATA2RAW
%
% Revision 1.23  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.22  2005/05/04 07:31:55  roboos
% remove avg after converting to single trial
%
% Revision 1.21  2005/03/08 10:39:50  roboos
% changed ordering: first check for backward compatibility, then set the (new) default
%
% Revision 1.20  2005/01/19 08:39:51  jansch
% added defaults for cfg.channel and cfg.channelcmb, delete cfg.channelcmb if output is not powandcsd
%
% Revision 1.19  2005/01/19 08:05:27  jansch
% fixed small bug in channelcombinations
%
% Revision 1.18  2005/01/18 15:20:38  roboos
% this function now already takes care of channelselection and, if needed, of selection of channelcombinations
%
% Revision 1.17  2005/01/18 15:05:15  roboos
% cleaned up configuration for sgn/label, now consistently using cfg.channel and cfg.channelcmb
%
% Revision 1.16  2004/11/01 11:34:53  roboos
% added support for timelocked trials, i.e. the result of timelockanalysis with keeptrials=yes
% tthese are now converted to raw trials prior to calling the freqanalysis_xxx subfunction
%
% Revision 1.15  2004/10/01 10:23:44  roboos
% fixed error in help for mtmconvol: f_timwin should be t_ftimwin
%
% Revision 1.14  2004/09/28 14:26:43  roboos
% fixed a typo in the help, added a line of comments
%
% Revision 1.13  2004/09/22 10:20:27  roboos
% converted to use external subfunctions time2offset and offset2time
% and add offset field to data structure if it is missing
%
% Revision 1.12  2004/09/21 12:07:28  marsie
% this version of the wrapper implements the new "freqanalysis_METHOD" convention
% for subfunction naming.
%
% Revision 1.11  2004/09/02 11:59:06  roboos
% restructured the help, fixed the copyrights
%
% Revision 1.10  2004/09/01 13:33:03  marsie
% switch to a wrapper function that calls the functions
% multitaperanalysis.m, wltanalysis.m and waveletanalysis.m for the corresponding
% methods
%

fieldtripdefs

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'raw', 'comp', 'mvar'}, 'feedback', 'yes', 'hasoffset', 'yes');

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'renamed',     {'label', 'channel'});
cfg = checkconfig(cfg, 'renamed',     {'sgn',   'channel'});
cfg = checkconfig(cfg, 'renamed',     {'labelcmb', 'channelcmb'});
cfg = checkconfig(cfg, 'renamed',     {'sgncmb',   'channelcmb'});
cfg = checkconfig(cfg, 'required',    {'method'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'fft',    'mtmfft'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'convol', 'mtmconvol'});

% select trials of interest
if ~isfield(cfg, 'trials'),   cfg.trials = 'all';  end % set the default
if ~strcmp(cfg.trials, 'all')
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
  data.offset = data.offset(cfg.trials);
  % update the trial definition (trl)
  if isfield(data, 'cfg')
    trl = findcfg(data.cfg, 'trl');
  else
    trl = [];
  end
  if isempty(trl)
    % a trial definition is expected in each continuous data set
    warning('could not locate the trial definition ''trl'' in the data structure');
  else
    cfg.trlold=trl;
    cfg.trl=trl(cfg.trials,:);
  end
end

% call the corresponding function
[freq] = feval(sprintf('freqanalysis_%s',lower(cfg.method)), cfg, data);
