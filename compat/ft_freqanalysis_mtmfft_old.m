function [freq] = ft_freqanalysis_mtmfft(cfg, data);

% FT_FREQANALYSIS_MTMFFT performs frequency analysis on any time series
% trial data using the 'multitaper method' (MTM) based on discrete
% prolate spheroidal sequences (Slepian sequences) as tapers. Alternatively,
% you can use conventional tapers (e.g. Hanning).
%
% Use as
%   [freq] = ft_freqanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration should be according to
%   cfg.method     = method used for frequency or time-frequency decomposition
%                    see FT_FREQANALYSIS for details
%   cfg.output     = 'pow'       return the power-spectra
%                    'powandcsd' return the power and the cross-spectra
%                    'fourier'   return the complex Fourier-spectra
%   cfg.taper      = 'dpss', 'hanning' or many others, see WINDOW (default = 'dpss')
%
% For cfg.output='powandcsd', you should specify the channel combinations
% between which to compute the cross-spectra as cfg.channelcmb. Otherwise
% you should specify only the channels in cfg.channel.
%
%   cfg.channel    = Nx1 cell-array with selection of channels (default = 'all'),
%                    see FT_CHANNELSELECTION for details
%   cfg.channelcmb = Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                    see FT_CHANNELCOMBINATION for details
%   cfg.foilim     = [begin end], frequency band of interest
%   cfg.tapsmofrq  = number, the amount of spectral smoothing through
%                    multi-tapering. Note that 4 Hz smoothing means
%                    plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrials = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.keeptapers = 'yes' or 'no', return individual tapers or average (default = 'no')
%   cfg.pad        = number or 'maxperlen', length in seconds to which the data can be padded out (default = 'maxperlen')
%
% The padding will determine your spectral resolution. If you want to
% compare spectra from data pieces of different lengths, you should use
% the same cfg.pad for both, in order to spectrally interpolate them to
% the same spectral resolution.  Note that this will run very slow if you
% specify cfg.pad as maxperlen AND the number of samples turns out to have
% a large prime factor sum. This is because the FFTs will then be computed
% very inefficiently.
%
% See also FT_FREQANALYSIS_MTMCONVOL, FT_FREQANALYSIS_WLTCONVOL, FT_FREQANALYSIS_TFR

% Undocumented local options
%   cfg.calcdof = 'yes'   calculate the degrees of freedom for every trial

% Copyright (c) 2003-2006, Pascal Fries, F.C. Donders Centre
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_freqanalysis_mtmfft.m 1974 2010-10-27 10:36:50Z jansch $

fieldtripdefs

% ensure that this function is started as a subfunction of the FT_FREQANALYSIS wrapper
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
  if ~strcmp(caller_name, 'ft_freqanalysis')
    error(['you should call FREQANALYSIS, instead of the ' upper(mfilename) ' subfunction']);
  end
end

% set all the defaults
if ~isfield(cfg, 'method'),        cfg.method     = 'mtmfft';               end
if ~isfield(cfg, 'keeptapers'),    cfg.keeptapers = 'no';                   end
if ~isfield(cfg, 'keeptrials'),    cfg.keeptrials = 'no';                   end
if ~isfield(cfg, 'calcdof'),       cfg.calcdof    = 'no';                   end
if ~isfield(cfg, 'pad'),           cfg.pad        = 'maxperlen';            end
if ~isfield(cfg, 'taper'),         cfg.taper      = 'dpss';                 end
if ~isfield(cfg, 'channel'),       cfg.channel    = 'all';                  end
if ~isfield(cfg, 'foilim'),        cfg.foilim     = [0 data.fsample/2];     end
if ~isfield(cfg, 'output'),
  if isfield(cfg, 'channelcmb') && ~isempty(cfg.channelcmb)
    cfg.output     = 'powandcsd';
  else
    cfg.output     = 'pow';
  end
end

if strcmp(cfg.output, 'fourier'),
  cfg.keeptrials = 'yes';
  cfg.keeptapers = 'yes';
end

if ~strcmp(cfg.method,'mtmfft')
  error('unsupported method');
end

% setting a flag (csdflg) that determines whether this routine outputs
% only power-spectra or power-spectra and cross-spectra?
if strcmp(cfg.output,'pow')
  powflg = 1;
  csdflg = 0;
  fftflg = 0;
elseif strcmp(cfg.output,'powandcsd')
  powflg = 1;
  csdflg = 1;
  fftflg = 0;
elseif strcmp(cfg.output,'fourier')
  powflg = 0;
  csdflg = 0;
  fftflg = 1;
else
  error('unsupported value for cfg.method');
end

if ~isfield(cfg, 'channelcmb') && csdflg
  %set the default for the channelcombination
  cfg.channelcmb = {'all' 'all'};
elseif isfield(cfg, 'channelcmb') && ~csdflg
  % no cross-spectrum needs to be computed, hence remove the combinations from cfg
  cfg = rmfield(cfg, 'channelcmb');
end

% ensure that channelselection and selection of channelcombinations is
% perfomed consistently
cfg.channel = ft_channelselection(cfg.channel, data.label);
if isfield(cfg, 'channelcmb')
  cfg.channelcmb = ft_channelcombination(cfg.channelcmb, data.label);
end

% determine the corresponding indices of all channels
sgnindx     = match_str(data.label, cfg.channel);
numsgn      = size(sgnindx,1);
if csdflg
  % determine the corresponding indices of all channel combinations
  sgncmbindx = zeros(size(cfg.channelcmb));
  for k=1:size(cfg.channelcmb,1)
    sgncmbindx(k,1) = find(strcmp(cfg.channelcmb(k,1), data.label));
    sgncmbindx(k,2) = find(strcmp(cfg.channelcmb(k,2), data.label));
    % this works the same, but is much slower in Octave
    % sgncmbindx(k,1) = strmatch(cfg.channelcmb(k,1), data.label, 'exact');
    % sgncmbindx(k,2) = strmatch(cfg.channelcmb(k,2), data.label, 'exact');
  end

  numsgncmb   = size(sgncmbindx,1);
  sgnindx     = unique([sgnindx(:); sgncmbindx(:)]);
  numsgn      = length(sgnindx);

  cutdatindcmb = zeros(size(sgncmbindx));
  for sgnlop = 1:numsgn
    cutdatindcmb(find(sgncmbindx == sgnindx(sgnlop))) = sgnlop;
  end
end

% if rectan is 1 it means that trials are of equal lengths
numper = numel(data.trial);
rectan = 1;
for perlop = 1:numper
  numdatbnsarr(perlop,1) = size(data.trial{perlop},2);
end
rectan = all(numdatbnsarr==numdatbnsarr(1));

% if cfg.pad is 'maxperlen', this is realized here:
if isequal(cfg.pad, 'maxperlen')
  cfg.pad = max(numdatbnsarr,[],1) ./ data.fsample;
else
  % check that the specified padding is not too short
  if cfg.pad<(max(numdatbnsarr,[],1)/data.fsample)
    error('the padding that you specified is shorter than the longest trial in the data');
  end
end
numsmp = ceil(cfg.pad .* data.fsample); % this used to be "cfg.pad .* data.fsample"

% keeping trials and/or tapers?
if strcmp(cfg.keeptrials,'no') &&  strcmp(cfg.keeptapers,'no')
  keeprpt = 1;
elseif strcmp(cfg.keeptrials,'yes') &&  strcmp(cfg.keeptapers,'no')
  keeprpt = 2;
elseif strcmp(cfg.keeptrials,'no') &&  strcmp(cfg.keeptapers,'yes')
  error('There is no support for keeping tapers WITHOUT KEEPING TRIALS.');
elseif strcmp(cfg.keeptrials,'yes') &&  strcmp(cfg.keeptapers,'yes')
  keeprpt = 4;
end

% calculating degrees of freedom
calcdof = strcmp(cfg.calcdof,'yes');

% doing the computation
boilim  = round(cfg.foilim ./ (data.fsample ./ numsmp)) + 1;
boi     = boilim(1):boilim(2);
numboi  = size(boi,2);
foi     = (boi-1) ./ cfg.pad;

if keeprpt == 1
  if powflg, powspctrm     = zeros(numsgn,numboi);              end
  if csdflg, crsspctrm     = complex(zeros(numsgncmb,numboi));  end
  if fftflg, fourierspctrm = complex(zeros(numsgn,numboi));     end
  dimord    = 'chan_freq';
elseif keeprpt == 2
  if powflg, powspctrm     = zeros(numper,numsgn,numboi);             end
  if csdflg, crsspctrm     = complex(zeros(numper,numsgncmb,numboi)); end
  if fftflg, fourierspctrm = complex(zeros(numper,numsgn,numboi));    end
  dimord    = 'rpt_chan_freq';
elseif keeprpt == 4
  if rectan == 1, % compute the amount of memory needed to collect the results
    numdatbns = numdatbnsarr(1,1);
    if strcmp(cfg.taper, 'dpss'),
      % ensure that the input arguments are double precision
      tap = double_dpss(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
    elseif strcmp(cfg.taper, 'sine')
      tap = sine_taper(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
    else
      tap(2,:) = nan;
    end
    numtap = size(tap,1)-1;
    numrpt = numtap.*numper;
  elseif rectan == 0,
    numrpt = 0;
    for perlop = 1:numper
      numdatbns = numdatbnsarr(perlop,1);
      if strcmp(cfg.taper, 'dpss'),
        % ensure that the input arguments are double precision
        tap = double_dpss(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
      elseif strcmp(cfg.taper, 'sine')
        tap = sine_taper(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
      else
        tap(2,:) = nan;
      end
      numtap = size(tap,1)-1;
      numrpt = numrpt + numtap;
    end
  end
  if powflg, powspctrm     = zeros(numrpt,numsgn,numboi);             end
  if csdflg, crsspctrm     = complex(zeros(numrpt,numsgncmb,numboi)); end
  if fftflg, fourierspctrm = complex(zeros(numrpt,numsgn,numboi));    end
  cnt = 0;
  dimord    = 'rpttap_chan_freq';
end

% these count the number of tapers
cumsumcnt = zeros(numper,1);
cumtapcnt = zeros(numper,1);

if rectan == 1
  % trials are of equal length, compute the set of tapers only once
  numdatbns = numdatbnsarr(1,1);
  if strcmp(cfg.taper, 'dpss')
    % create a sequence of DPSS (Slepian) tapers
    % ensure that the input arguments are double precision
    tap = double_dpss(numdatbns,numdatbns*(cfg.tapsmofrq./data.fsample))';
  elseif strcmp(cfg.taper, 'sine')
    tap = sine_taper(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
  else
    % create a single taper according to the window specification as a
    % replacement for the DPSS (Slepian) sequence
    tap = window(cfg.taper, numdatbns)';
    tap = tap./norm(tap);
    % freqanalysis_mtmfft always throws away the last taper of the Slepian sequence, so add a dummy taper
    tap(2,:) = nan;
  end
  numtap = size(tap,1) - 1;
  if keeprpt == 2 || calcdof
    cumtapcnt(:) = numtap;
  end
  if (numtap < 1)
    error(sprintf(...
      'datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',...
      numdatbns/data.fsample, cfg.tapsmofrq, data.fsample/numdatbns));
  elseif (numtap < 2) && strcmp(cfg.taper, 'dpss')
    fprintf('WARNING: using only one taper for specified smoothing\n');
  end
  pad = zeros(1,numsmp - numdatbns);
end

for perlop = 1:numper
  fprintf('processing trial %d, ', perlop);
  if keeprpt == 2
    cnt = perlop;
    cumsumcnt(cnt,1) = numdatbnsarr(perlop,1);
  end
  if rectan == 0
    % trials are not of equal length, compute the set of tapers for this trial
    numdatbns = numdatbnsarr(perlop,1);
    if strcmp(cfg.taper, 'dpss')
      % create a sequence of DPSS (Slepian) tapers
      % ensure that the input arguments are double precision
      tap = double_dpss(numdatbns,numdatbns*(cfg.tapsmofrq./data.fsample))';
    elseif strcmp(cfg.taper, 'sine')
      tap = sine_taper(numdatbns, numdatbns*(cfg.tapsmofrq./data.fsample))';
    else
      % create a single taper according to the window specification as a
      % replacement for the DPSS (Slepian) sequence
      tap = window(cfg.taper, numdatbns)';
      tap = tap./norm(tap);
      % freqanalysis_mtmfft always throws away the last taper of the Slepian sequence, so add a dummy taper
      tap(2,:) = nan;
    end
    numtap = size(tap,1) - 1;
    if keeprpt == 2 || calcdof
      cnt = perlop;
      cumtapcnt(cnt,1) = numtap;
    end
    if (numtap < 1)
      error(sprintf(...
        'datalength to short for specified smoothing\ndatalength: %.3f s, smoothing: %.3f Hz, minimum smoothing: %.3f Hz',...
        numdatbns/data.fsample, cfg.tapsmofrq, data.fsample/numdatbns));
    elseif (numtap < 2) && strcmp(cfg.taper, 'dpss')
      fprintf('WARNING: using only one taper for specified smoothing\n');
    end
    pad = zeros(1,numsmp - numdatbns);
  end
  for taplop = 1:numtap
    if keeprpt == 4
      cnt = cnt+1;
      cumsumcnt(perlop,1) = numdatbnsarr(perlop,1);
      cumtapcnt(perlop,1) = numtap;
    end
    if calcdof
      cumtapcnt(perlop,1) = numtap;
    end

    autspctrmacttap = complex(zeros(numsgn,numboi));
    for sgnlop = 1:numsgn
      dum = fft([data.trial{perlop}(sgnindx(sgnlop),:) .* tap(taplop,:) , ...
        pad],[],2);
      autspctrmacttap(sgnlop,:) = dum(boi);
    end
    if taplop == 1
      fprintf('nfft: %d samples, taper length: %d samples, %d tapers\n',length(dum),size(tap,2),numtap);
    end
    if powflg
      powdum = 2 .* (autspctrmacttap .* conj(autspctrmacttap)) ./ numsmp; %cf Numercial Receipes 13.4.9
      if keeprpt == 1
        powspctrm(:,:) = powspctrm(:,:) + (powdum ./ numtap);
      elseif keeprpt == 2
        powspctrm(cnt,:,:) = powspctrm(cnt,:,:) + (permute(powdum,[3,1,2]) ./ numtap);
      elseif keeprpt == 4
        powspctrm(cnt,:,:) = powdum;
      end
    end
    if fftflg
      fourierdum = (autspctrmacttap) .* sqrt(2 ./ numsmp); %cf Numercial Receipes 13.4.9
      if keeprpt == 1
        fourierspctrm(:,:) = fourierspctrm(:,:) + (fourierdum ./ numtap);
      elseif keeprpt == 2
        fourierspctrm(cnt,:,:) = fourierspctrm(cnt,:,:) + (permute(fourierdum,[3,1,2]) ./ numtap);
      elseif keeprpt == 4
        fourierspctrm(cnt,:,:) = fourierdum;
      end
    end
    if csdflg
      csddum = 2.* (autspctrmacttap(cutdatindcmb(:,1),:) .* ...
        conj(autspctrmacttap(cutdatindcmb(:,2),:))) ./ numsmp;
      if keeprpt == 1
        crsspctrm(:,:) = crsspctrm(:,:) + csddum ./ numtap;
      elseif keeprpt == 2
        crsspctrm(cnt,:,:) = crsspctrm(cnt,:,:) + permute(csddum,[3,1,2]) ./ numtap;
      elseif keeprpt == 4
        crsspctrm(cnt,:,:) = csddum;
      end
    end
  end % taplop
end % perlop
if keeprpt ==1
  if powflg, powspctrm = powspctrm ./ numper; end
  if csdflg, crsspctrm = crsspctrm ./ numper; end
end

% collect the results
freq.label      = data.label(sgnindx);
freq.dimord     = dimord;
freq.freq       = foi;
hasdc           = find(freq.freq==0); 

if powflg
  % correct the 0 Hz bin if present, scaling with a factor of 2 is only appropriate for ~0 Hz
  if ~isempty(hasdc)
    if keeprpt>1
      powspctrm(:,:,hasdc,:) = powspctrm(:,:,hasdc,:)./2;      
    else
      powspctrm(:,hasdc,:) = powspctrm(:,hasdc,:)./2;
    end
  end
  freq.powspctrm  = powspctrm;
end

if fftflg
  % correct the 0 Hz bin if present
  if ~isempty(hasdc)
    if keeprpt>1
      fourierspctrm(:,:,hasdc,:) = fourierspctrm(:,:,hasdc,:)./sqrt(2);      
    else
      fourierspctrm(:,hasdc,:) = fourierspctrm(:,hasdc,:)./sqrt(2);
    end
  end
  freq.fourierspctrm  = fourierspctrm;
end

if csdflg
  % correct the 0 Hz bin if present
  if ~isempty(hasdc)
    if keeprpt>1
      crsspctrm(:,:,hasdc,:) = crsspctrm(:,:,hasdc,:)./2;      
    else
      crsspctrm(:,hasdc,:) = crsspctrm(:,hasdc,:)./2;
    end
  end
  freq.labelcmb   = cfg.channelcmb;
  freq.crsspctrm  = crsspctrm;
end

if strcmp(cfg.method,'mtmfft') && (keeprpt == 2 || keeprpt == 4)
  freq.cumsumcnt = cumsumcnt;
end

if strcmp(cfg.method,'mtmfft') && (keeprpt == 2 || keeprpt == 4)
  freq.cumtapcnt = cumtapcnt;
end

if calcdof
  freq.dof=2*repmat(cumtapcnt,[1,numboi]);
end;

try, freq.grad = data.grad; end   % remember the gradiometer array
try, freq.elec = data.elec; end   % remember the electrode array

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add information about the version of this function to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i1] = dbstack;
  cfg.version.name = st(i1);
end
cfg.version.id = '$Id: ft_freqanalysis_mtmfft.m 1974 2010-10-27 10:36:50Z jansch $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for Matlab 6.5 and 7.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tap] = double_dpss(a, b, varargin);
tap = dpss(double(a), double(b), varargin{:});

