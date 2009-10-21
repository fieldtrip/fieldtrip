function [freq] = freqanalysis_tfr(cfg, data);

% FREQANALYSIS_TFR computes time-frequency representations of single-trial
% data using a convolution in the time-domain with Morlet's wavelets.
%
% Use as
%   [freq] = freqanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from
% the PREPROCESSING function. The configuration should be according to
%   cfg.method        = method used for frequency or time-frequency decomposition
%                       see FREQANALYSIS for details
%   cfg.foi           = vector 1 x numfoi, frequencies of interest
%   cfg.waveletwidth  = 'width' of wavelets expressed in cycles (default = 7)
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see CHANNELSELECTION for details
%   cfg.downsample    = ratio for downsampling, which occurs after convolution (default = 1)
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrials    = 'yes' or 'no', return individual trials or average (default = 'no')
%
% See also FREQANALYSIS

% Undocumented local options:
% cfg.latency
% cfg.output

% Copyright (C) 2003, Ole Jensen, FCDC
%
% $Log: freqanalysis_tfr.m,v $
% Revision 1.23  2009/02/16 20:30:44  jansch
% changed inconsistent normalization of power. initially, normalization was
% (2./data.fsample).^2. in other freqanalysis_xxx functions, this normalization
% is (2./data.fsample)
%
% Revision 1.22  2008/11/11 18:59:26  sashae
% added call to checkconfig at end of function (trackconfig and checksize)
%
% Revision 1.21  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.20  2008/01/18 13:14:50  sashae
% added option for trial selection, updated documentation
%
% Revision 1.19  2007/03/27 11:00:28  ingnie
% deleted call to fixdimord because this is low level function
%
% Revision 1.18  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.17  2006/06/06 16:55:45  ingnie
% updated documentation, deleted default cfg.foi
%
% Revision 1.16  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.15  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.14  2006/03/29 08:44:28  roboos
% added support for keeptrials, added check for call from
% freqanalysis-wrapper, changed default feedback, cleaned up help,
% autoindented code (spaces/tabs).
%
% Revision 1.13  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.12  2006/02/07 20:08:22  roboos
% changed all occurences of a dimord with chancmb (was previous sgncmb) into chan
%
% Revision 1.11  2006/02/07 08:12:19  roboos
% changed fprintf into progress indicator, added cfg.feedback
%
% Revision 1.10  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.9  2005/08/05 09:14:03  roboos
% replaced computation of trialtime by the time axis that is present in the data
% added warning to use RESAMPLEDATA is cfg.downsample is used
%
% Revision 1.8  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.7  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.6  2005/01/19 08:40:40  jansch
% added check and error message if output=powandcsd
%
% Revision 1.5  2005/01/18 15:11:39  roboos
% Cleaned up configuration for sgn/sgncmb, now exclusively using channel/channelcmb which is consistent with rest of fieldtrip and freqanalysis documentation. Also the output now only contains freq.label/labelcmb and not any more the old sgn/sgncmb.
%
% Revision 1.4  2004/11/01 11:33:37  roboos
% replaced obbsolete translate_channel_list by channelselection
%
% Revision 1.3  2004/10/28 09:00:33  roboos
% fixed bug in caller function detection (for matlab 6.1 and 6.5)
%
% Revision 1.2  2004/10/28 07:21:46  roboos
% added check to ensure that the FREQANALYSIS wrapper is started instead of the
% FERQANALYSIS_xxx subfunctions
%
% Revision 1.1  2004/09/21 12:05:26  marsie
% the tfr method of waveletanalysis.m has been moved to this seperate function
%
% Revision 1.5  2004/08/16 13:31:58  roboos
% added a ";" at the end of the line
%
% Revision 1.4  2004/07/01 11:57:44  olejen
% Error in dimord corrected
%
% Revision 1.3  2004/06/30 11:44:48  roboos
% added CVS Log feature and copyright statement, no code changes
% (note: the previous change from revision 1.1 to 1.2 only affected
% documentation and comments, and also did not include any change to
% the functionality of the code)

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

% set the defaults
if ~isfield(cfg, 'method'),         cfg.method  = 'tfr';          end
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';          end
if ~isfield(cfg, 'latency'),        cfg.latency = 'minperlength'; end
if ~isfield(cfg, 'keeptrials'),     cfg.keeptrials   = 'no';      end
if ~isfield(cfg, 'waveletwidth'),   cfg.waveletwidth = 7;         end
if ~isfield(cfg, 'downsample'),     cfg.downsample   = 1;         end
if ~isfield(cfg, 'feedback'),       cfg.feedback     = 'text';    end

if isfield(cfg, 'output') && strcmp(cfg.output, 'powandcsd'),
  error('This function does not compute cross-spectra\n');
end

% determine the channels of interest
cfg.channel = channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);

% determine the duration of each trial
ntrial = length(data.trial);
nchan = size(data.trial{1}, 1);

for i=1:ntrial
  nsampl(i)         = size(data.trial{i}, 2);
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
end;

if cfg.downsample > 1
  % perform a decimation of the input data
  warning('decimating the input data, better is to use RESAMPLEDATA');
  for k=1:ntrial
    dTmp = data.trial{k};
    data.trial{k} = dTmp(:,1:cfg.downsample:end);
    tTmp = data.time{k};
    data.time{k} = tTmp(1:cfg.downsample:end);
  end
  data.fsample = data.fsample / cfg.downsample;
end

% automatically determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];

% latency window for averaging and variance computation is given in seconds
if (strcmp(cfg.latency, 'minperlength'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = minperlength(2);
end
if (strcmp(cfg.latency, 'prestim'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = 0;
end
if (strcmp(cfg.latency, 'poststim'))
  cfg.latency = [];
  cfg.latency(1) = 0;
  cfg.latency(2) = minperlength(2);
end

M = waveletfam(cfg.foi,data.fsample,cfg.waveletwidth);

progress('init', cfg.feedback, 'convolving wavelets');

for i=1:ntrial
  indicvect = data.time{i};
  progress(i/ntrial, 'convolving wavelets, trial %d of %d\n', i, ntrial);

  %for average and variance
  begsampl = nearest(indicvect,cfg.latency(1));
  endsampl = nearest(indicvect,cfg.latency(2));

  numsamples(i) = endsampl-begsampl+1;

  if (i==1)
    % allocate memory to hold the resulting powerspectra
    if strcmp(cfg.keeptrials, 'yes')
      freq.powspctrm = zeros(ntrial,nchan,length(cfg.foi),ceil((endsampl-begsampl+1)/cfg.downsample));
    else
      freq.powspctrm = zeros(nchan,length(cfg.foi),ceil((endsampl-begsampl+1)/cfg.downsample));
    end
  end;

  dat = data.trial{i}(chansel,begsampl:endsampl);
  for k=1:size(dat,1)
    for j=1:length(cfg.foi)
      cTmp = conv(dat(k,:),M{j});
      cTmp = 2*(abs(cTmp).^2)/data.fsample;
      cTmp = cTmp(ceil(length(M{j})/2):length(cTmp)-floor(length(M{j})/2));
      cTmp = cTmp(:,1:cfg.downsample:end);
      if strcmp(cfg.keeptrials, 'yes')
        freq.powspctrm(i,k,j,:) = cTmp';
      else
        freq.powspctrm(k,j,:) = squeeze(freq.powspctrm(k,j,:)) + cTmp';  % compute the running sum
      end
    end
  end

end %for ntrial

progress('close');

if strcmp(cfg.keeptrials, 'yes')
  freq.dimord    = 'rpt_chan_freq_time';
else
  freq.dimord    = 'chan_freq_time';
  freq.powspctrm = freq.powspctrm / ntrial;  % compute the average
end
freq.label     = cfg.channel;
freq.freq      = cfg.foi;
freq.time      = indicvect(1:cfg.downsample:end);

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: freqanalysis_tfr.m,v 1.23 2009/02/16 20:30:44 jansch Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for waveletanalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = waveletfam(foi,Fs,width)
dt = 1/Fs;
for k=1:length(foi)
  sf = foi(k)/width;
  st = 1/(2*pi*sf);
  toi=-3.5*st:dt:3.5*st;
  A = 1/sqrt(st*sqrt(pi));
  M{k}= A*exp(-toi.^2/(2*st^2)).*exp(i*2*pi*foi(k).*toi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for waveletanalysis
%
% Return a vector containing the energy as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets.
% s : signal
% Fs: sampling frequency
% width: width of Morlet wavelet (>= 5 suggested).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = energyvec(f,s,Fs,width)
dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);
t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
size(m)
y = conv(s,m);
y = (2*abs(y)/Fs).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction for waveletanalysis
%
% Morlet's wavelet for frequency f and time t.
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet.
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY
%
% Ole Jensen, August 1998
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = morlet(f,t,width)
sf = f/width;
st = 1/(2*pi*sf);
A = 1/sqrt(st*sqrt(pi));
y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);
