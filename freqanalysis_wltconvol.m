function [freq] = freqanalysis_wltconvol(cfg, data);

% FREQANALYSIS_WLTCONVOL performs time-frequency analysis on any time series trial data
% using the 'wavelet method' based on Morlet wavelets.
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
%
% For cfg.output='powandcsd', you should specify the channel combinations
% between which to compute the cross-spectra as cfg.channelcmb. Otherwise
% you should specify only the channels in cfg.channel.
%
%   cfg.channel    = Nx1 cell-array with selection of channels (default = 'all'),
%                    see CHANNELSELECTION for details
%   cfg.channelcmb = Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                    see CHANNELCOMBINATION for details
%   cfg.foi        = vector 1 x numfoi, frequencies of interest
%   cfg.toi        = vector 1 x numtoi, the times on which the analysis windows
%                    should be centered (in seconds)
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrials = 'yes' or 'no', return individual trials or average (default = 'no')
%   cfg.width      = 'width' of the wavelet, determines the temporal and spectral
%                    resolution of the analysis (default = 7)
%                    constant, for a 'classical constant-Q' wavelet analysis
%                    vector, defining a variable width for each frequency
%   cfg.gwidth     = determines the length of the used wavelets in standard deviations
%                    of the implicit Gaussian kernel and should be choosen
%                    >= 3; (default = 3)
%
% The standard deviation in the frequency domain (sf) at frequency f0 is
% defined as: sf = f0/width
% The standard deviation in the temporal domain (st) at frequency f0 is
% defined as: st = width/f0 = 1/sf
%
% See also FREQANALYSIS

% Copyright (C) 2003-2007, Markus Siegel, F.C. Donders Centre
%
% $Log: freqanalysis_wltconvol.m,v $
% Revision 1.22  2008/11/11 18:59:26  sashae
% added call to checkconfig at end of function (trackconfig and checksize)
%
% Revision 1.21  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.20  2008/05/13 15:36:09  roboos
% fixed potential bug in assessing the number of trials (when data.trial was column instead of row vector)a, now use numel instead of size
%
% Revision 1.19  2008/01/18 13:14:50  sashae
% added option for trial selection, updated documentation
%
% Revision 1.18  2007/08/06 15:00:52  roboos
% updated documentaton and copyright
%
% Revision 1.17  2007/06/14 12:43:12  jansch
% fixed improper normalisation of powdum
%
% Revision 1.16  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.15  2006/06/20 16:28:29  ingnie
% added consistent handling of cfg.channel and cfg.channelcmb
%
% Revision 1.14  2006/06/13 14:53:22  ingnie
% some change in white space defaults, added default cfg.channel = 'all'
%
% Revision 1.13  2006/06/06 16:57:51  ingnie
% updated documentation
%
% Revision 1.12  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.11  2006/02/07 20:08:22  roboos
% changed all occurences of a dimord with chancmb (was previous sgncmb) into chan
%
% Revision 1.10  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.9  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.8  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.7  2005/01/19 08:42:42  jansch
% removed obsolete code for generating and checking channelcombinations
% cleaned up handling of sgnindx/sgncmbindx
% ensured that power is computed for all channels that are in channelcmb
%
% Revision 1.6  2005/01/18 15:11:39  roboos
% Cleaned up configuration for sgn/sgncmb, now exclusively using channel/channelcmb which is consistent with rest of fieldtrip and freqanalysis documentation. Also the output now only contains freq.label/labelcmb and not any more the old sgn/sgncmb.
%
% Revision 1.5  2004/12/20 14:52:56  roboos
% changed rounding off of timboi
%
% Revision 1.4  2004/11/18 16:51:58  roboos
% updated the help: pointed to CHANNELCOMBINATION for cfg.channelcmb
%
% Revision 1.3  2004/10/28 09:00:33  roboos
% fixed bug in caller function detection (for matlab 6.1 and 6.5)
%
% Revision 1.2  2004/10/28 07:21:46  roboos
% added check to ensure that the FREQANALYSIS wrapper is started instead of the
% FERQANALYSIS_xxx subfunctions
%
% Revision 1.1  2004/09/21 12:04:53  marsie
% the wltconvol method of wltanalysis.m has been moved to this seperate function
%
% Revision 1.2  2004/08/25 20:47:18  marsie
% fixed scaling problems
%
% Revision 1.1  2004/08/25 19:18:43  marsie
% initial release
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

% set all the defaults
if ~isfield(cfg, 'method'),        cfg.method     = 'wltconvol';  end
if ~isfield(cfg, 'keeptrials'),    cfg.keeptrials = 'no';         end
if ~isfield(cfg, 'output'),        cfg.output     = 'powandcsd';  end
if ~isfield(cfg, 'pad'),           cfg.pad        = 'maxperlen';  end
if ~isfield(cfg, 'width'),         cfg.width      = 7;            end
if ~isfield(cfg, 'gwidth'),        cfg.gwidth     = 3;            end
if ~isfield(cfg, 'channel'),       cfg.channel    = 'all';        end

% expand cfg.width to array if constant width
if prod(size(cfg.width)) == 1
  cfg.width = ones(1,length(cfg.foi)) * cfg.width;
end

% setting a flag (csdflg) that determines whether this routine outputs
% only power-spectra or power-spectra and cross-spectra?
if strcmp(cfg.output,'pow')
  csdflg = 0;
elseif strcmp(cfg.output,'powandcsd')
  csdflg = 1;
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
cfg.channel = channelselection(cfg.channel, data.label);
if isfield(cfg, 'channelcmb')
  cfg.channelcmb = channelcombination(cfg.channelcmb, data.label);
end

% determine the corresponding indices of all channels
sgnindx     = match_str(data.label, cfg.channel);
numsgn      = size(sgnindx,1);
if csdflg
  % determine the corresponding indices of all channel combinations
  for k=1:size(cfg.channelcmb,1)
    sgncmbindx(k,1) = strmatch(cfg.channelcmb(k,1), data.label, 'exact');
    sgncmbindx(k,2) = strmatch(cfg.channelcmb(k,2), data.label, 'exact');
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
  if numdatbnsarr(perlop,1) ~= numdatbnsarr(1,1)
    rectan = 0;
  end
end

%if cfg.pad is 'maxperlen', this is realized here:
if ischar(cfg.pad)
  if strcmp(cfg.pad,'maxperlen')
    cfg.pad = max(numdatbnsarr,[],1) ./ data.fsample;
  end
end
numsmp = cfg.pad .* data.fsample;

% keeping trials and/or tapers?
if strcmp(cfg.keeptrials,'no')
  keep = 1;
elseif strcmp(cfg.keeptrials,'yes')
  keep = 2;
end

% do the computation for WLTCONVOL
if strcmp(cfg.method,'wltconvol')
  minoffset = min(data.offset);
  timboi = round(cfg.toi .* data.fsample - minoffset);
  toi = round(cfg.toi .* data.fsample) ./ data.fsample;
  numtoi = length(toi);
  numfoi = length(cfg.foi);
  knlspctrmstr = cell(numfoi,1);
  for foilop = 1:numfoi
    dt = 1/data.fsample;
    sf = cfg.foi(foilop)/cfg.width(foilop);
    st = 1/(2*pi*sf);
    toi2 = -cfg.gwidth*st:dt:cfg.gwidth*st;
    A = 1/sqrt(st*sqrt(pi));
    tap = (A*exp(-toi2.^2/(2*st^2)))';
    acttapnumsmp = size(tap,1);
    taplen(foilop) = acttapnumsmp;
    ins = ceil(numsmp./2) - floor(acttapnumsmp./2);
    prezer = zeros(ins,1);
    pstzer = zeros(numsmp - ((ins-1) + acttapnumsmp)-1,1);
    ind = (0:acttapnumsmp-1)' .* ...
      ((2.*pi./data.fsample) .* cfg.foi(foilop));
    knlspctrmstr{foilop} = complex(zeros(1,numsmp));
    knlspctrmstr{foilop} = ...
      fft(complex(vertcat(prezer,tap.*cos(ind),pstzer), ...
      vertcat(prezer,tap.*sin(ind),pstzer)),[],1)';
  end
  if keep == 1
    powspctrm = zeros(numsgn,numfoi,numtoi);
    if csdflg, crsspctrm = complex(zeros(numsgncmb,numfoi,numtoi)); end
    cntpertoi = zeros(numfoi,numtoi);
    dimord    = 'chan_freq_time';
  elseif keep == 2
    powspctrm = zeros(numper,numsgn,numfoi,numtoi);
    if csdflg, crsspctrm = complex(zeros(numper,numsgncmb,numfoi,numtoi)); end
    dimord    = 'rpt_chan_freq_time';
  end
  for perlop = 1:numper
    fprintf('processing trial %d: %d samples\n', perlop, numdatbnsarr(perlop,1));
    if keep == 2
      cnt = perlop;
    end
    numdatbns = numdatbnsarr(perlop,1);
    prepad = zeros(1,data.offset(perlop) - minoffset);
    pstpad = zeros(1,minoffset + numsmp - (data.offset(perlop) + numdatbns));
    datspctra = complex(zeros(numsgn,numsmp));
    for sgnlop = 1:numsgn
      datspctra(sgnlop,:) = fft([prepad, data.trial{perlop}(sgnindx(sgnlop),:), ...
        pstpad],[],2);
    end
    for foilop = 1:numfoi
      fprintf('processing frequency %d (%.2f Hz)\n', foilop,cfg.foi(foilop));
      actfoinumsmp = taplen(foilop);
      acttimboiind = ...
        find(timboi >= (-minoffset + data.offset(perlop) + (actfoinumsmp ./ 2)) & ...
        timboi < (-minoffset + data.offset(perlop) + numdatbns - (actfoinumsmp ./2)));
      nonacttimboiind = ...
        find(timboi < (-minoffset + data.offset(perlop) + (actfoinumsmp ./ 2)) | ...
        timboi >= (-minoffset + data.offset(perlop) + numdatbns - (actfoinumsmp ./2)));
      acttimboi = timboi(acttimboiind);
      numacttimboi = length(acttimboi);
      if keep ==1
        cntpertoi(foilop,acttimboiind) = cntpertoi(foilop,acttimboiind) + 1;
      end
      if keep == 3
        cnt = 1;
      end
      autspctrmacttap = complex(zeros(numsgn,numacttimboi));
      if numacttimboi > 0
        for sgnlop = 1:numsgn
          dum = fftshift(ifft(datspctra(sgnlop,:) .* ...
            [knlspctrmstr{foilop}],[],2));
          autspctrmacttap(sgnlop,:) = dum(acttimboi);
        end
        powdum = 2.* (abs(autspctrmacttap).^2)  ./ data.fsample;
      end
      if keep == 1 && numacttimboi > 0
        powspctrm(:,foilop,acttimboiind) = powspctrm(:,foilop,acttimboiind) + ...
          reshape(powdum,[numsgn,1,numacttimboi]);
      elseif keep == 2 && numacttimboi > 0
        powspctrm(cnt,:,foilop,acttimboiind) = powspctrm(cnt,:,foilop,acttimboiind) + ...
          reshape(powdum,[1,numsgn,1,numacttimboi]);
        powspctrm(cnt,:,foilop,nonacttimboiind) = nan;
      elseif keep == 2 && numacttimboi == 0
        powspctrm(cnt,:,foilop,nonacttimboiind) = nan;
      end
      if csdflg
        csddum = 2.* (autspctrmacttap(cutdatindcmb(:,1),:) .* ...
          conj(autspctrmacttap(cutdatindcmb(:,2),:))) ./ data.fsample; %actfoinumsmp;
        if keep == 1 && numacttimboi > 0
          crsspctrm(:,foilop,acttimboiind) = ...
            crsspctrm(:,foilop,acttimboiind) + ...
            reshape(csddum,[numsgncmb,1,numacttimboi]);
        elseif keep == 2 && numacttimboi > 0
          crsspctrm(cnt,:,foilop,acttimboiind) = ...
            crsspctrm(cnt,:,foilop,acttimboiind) + ...
            reshape(csddum,[1,numsgncmb,1,numacttimboi]);
          crsspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        elseif keep == 2 && numacttimboi == 0
          crsspctrm(cnt,:,foilop,nonacttimboiind) = nan;
        end
      end
    end% of foilop
  end%of perlop
  if keep == 1
    warning off
    powspctrm(:,:,:) = powspctrm(:,:,:) ./ repmat(permute(cntpertoi,[3,1,2]),[numsgn,1,1]);
    if csdflg
      crsspctrm(:,:,:) = crsspctrm(:,:,:) ./ repmat(permute(cntpertoi,[3,1,2]),[numsgncmb,1,1]);
    end
    warning on
  end
end

% collect the results
freq.label      = data.label(sgnindx);
freq.dimord     = dimord;
freq.powspctrm  = powspctrm;
freq.freq       = cfg.foi;
freq.time       = cfg.toi;

if csdflg
  freq.labelcmb   = cfg.channelcmb;
  freq.crsspctrm  = crsspctrm;
end

try, freq.grad = data.grad; end   % remember the gradiometer array
try, freq.elec = data.elec; end   % remember the electrode array

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add information about the version of this function to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i1] = dbstack;
  cfg.version.name = st(i1);
end
cfg.version.id = '$Id: freqanalysis_wltconvol.m,v 1.22 2008/11/11 18:59:26 sashae Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

