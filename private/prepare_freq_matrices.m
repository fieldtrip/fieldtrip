function [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices(cfg, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that converts a freq structure into Cf, Cr and Pr
% this is used in sourecanalysis
%
% This function returns data matrices with a channel order that is consistent
% with the original channel order in the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2004-2006, Robert Oostenveld
%
% $Log: prepare_freq_matrices.m,v $
% Revision 1.25  2009/08/16 12:43:22  jansch
% added default empty cfg.refchan
%
% Revision 1.24  2009/02/05 10:22:07  roboos
% better support for single-trial data, thanks to Vladimir
%
% Revision 1.23  2008/07/11 14:06:43  jansch
% fixed the cfg.channel freq.label issue for fourierspctrm as well
%
% Revision 1.22  2008/07/10 15:50:27  roboos
% keep channel order consistent with the input cfg
%
% Revision 1.21  2008/07/08 15:39:22  roboos
% initial version for Saskia to work on
%
% Revision 1.20  2006/10/30 16:37:29  roboos
% fixed bug in reshaping of fourier data (was ctranspose, whereas it should be normal transpose), all phase differences in the CSD matix were therefore defined the other way around
% added support for Cr and Pr (needed for dics) in case fourier input data
%
% Revision 1.19  2006/10/02 13:02:27  roboos
% fixed bug for rpttap_chan_freq_time, where the 'tap' part was dropped from the dimord
%
% Revision 1.18  2006/03/08 15:31:14  roboos
% removed default cfg.channel=all, since no channelselection is done inside this function
%
% Revision 1.17  2006/02/23 10:28:17  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.16  2006/02/10 09:15:30  jansch
% changed some chancmb which were missed during last fix
%
% Revision 1.15  2006/02/10 09:12:15  jansch
% changed _chancmb in the dimord into _chan
%
% Revision 1.14  2006/02/01 12:26:04  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.13  2005/08/26 14:56:27  roboos
% selection of tbin in freq.toi was done twice, resulting in error
%
% Revision 1.12  2005/08/26 13:39:02  roboos
% changed the indentation throughout the code, no functional changes
%
% Revision 1.11  2005/08/25 08:05:42  jansch
% removed support for fourier in case of average, since that does not make sense
% added support for multiple tapers in computation of CSD from fourier input
%
% Revision 1.10  2005/08/16 13:16:24  jansch
% *** empty log message ***
%
% Revision 1.9  2005/08/16 12:47:57  jansch
% made change to also support input dimord of 'rpttap_sgncmb_frq'
%
% Revision 1.8  2005/08/15 13:30:50  jansch
% removed the normalisation for the number of samples, when computing the cross-
% spectra from the fourierspctrum: TAKE CARE! this is not backward-compatible
%
% Revision 1.7  2005/05/20 16:50:26  roboos
% added normalisation of the CSD matrix computed from the fourier spectra
%
% Revision 1.6  2005/05/17 18:01:21  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
% added cfg.dicsfix for taking conjugate of csd from MEG to reference channel (default)
% implemented computation of csd from single trial fourier spectrum (if freqanalysis with output=fourier)
%
% Revision 1.5  2005/03/03 10:59:29  roboos
% Added cfg as output argument and added selected chanels to the output
% cfg, so that the calling function knows which were selected. Renamed
% tbin into timebin.
%
% Revision 1.4  2005/02/16 15:13:50  roboos
% renamed cfg.refchannel into cfg.refchan, replaced detection of coh_refchan submethod for dics (now looking for presence of cfg.refchan instead of explicitely testing the submethod)
%
% Revision 1.3  2004/08/26 12:43:22  roboos
% fixed bug in backward-compatibility code for sgn/sgncmb
%
% Revision 1.2  2004/08/25 17:00:08  roboos
% changed from sgn/sgncmb to label/labelcmb, backward compatible
%
% Revision 1.1  2004/06/28 08:59:38  roboos
% moved files from fieldtrip to fieldtrip/private
%
% Revision 1.1  2004/04/08 15:51:20  roberto
% initial submissiion into cvs
%

% set the defaults
if ~isfield(cfg, 'dicsfix'), cfg.dicsfix = 'yes'; end
if ~isfield(cfg, 'quickflag'), cfg.quickflag = 0; end
if ~isfield(cfg, 'refchan'), cfg.refchan = []; end

quickflag = cfg.quickflag==1;

Cf = [];
Cr = [];
Pr = [];

% the old freqanalysis (up to revision 1.9) used sgn and sgncmb
% for backward compatibility rename these
if isfield(freq, 'sgn')
  freq.label = freq.sgn;
  freq = rmfield(freq, 'sgn');
end
if isfield(freq, 'sgncmb')
  freq.labelcmb = freq.sgncmb;
  freq = rmfield(freq, 'sgncmb');
end

% select the latency of interest for time-frequency data
if strcmp(freq.dimord, 'chan_freq_time')
  tbin = nearest(freq.time, cfg.latency);
  fprintf('selecting timeslice %d\n', tbin);
  freq.time = freq.time(tbin);
  % remove all other latencies from the data structure and reduce the number of dimensions
  if isfield(freq, 'powspctrm'),     freq.powspctrm     = squeeze(freq.powspctrm(:,:,tbin));     end;
  if isfield(freq, 'crsspctrm'),     freq.crsspctrm     = squeeze(freq.crsspctrm(:,:,tbin));     end;
  if isfield(freq, 'fourierspctrm'), freq.fourierspctrm = squeeze(freq.fourierspctrm(:,:,tbin)); end;
  freq.dimord = freq.dimord(1:(end-5));  % remove the '_time' part
elseif strcmp(freq.dimord, 'rpt_chan_freq_time') || strcmp(freq.dimord, 'rpttap_chan_freq_time')
  tbin = nearest(freq.time, cfg.latency);
  fprintf('selecting timeslice %d\n', tbin);
  freq.time = freq.time(tbin);
  % remove all other latencies from the data structure and reduce the number of dimensions
  if isfield(freq, 'powspctrm'),    freq.powspctrm     = squeeze(freq.powspctrm(:,:,:,tbin));      end;
  if isfield(freq, 'crsspctrm'),    freq.crsspctrm     = squeeze(freq.crsspctrm(:,:,:,tbin));      end;
  if isfield(freq, 'fourierspctrm') freq.fourierspctrm = squeeze(freq.fourierspctrm(:,:,:,tbin));  end;
  freq.dimord = freq.dimord(1:(end-5));  % remove the '_time' part
else
  tbin = [];
end

% the time-frequency latency has already been squeezed away (see above)
if strcmp(freq.dimord, 'chan_freq')
  Ntrials = 1;
elseif strcmp(freq.dimord, 'rpt_chan_freq')
  Ntrials = length(freq.cumtapcnt);
elseif strcmp(freq.dimord, 'rpttap_chan_freq')
  Ntrials = length(freq.cumtapcnt);
else
  error('unrecognized dimord for frequency data');
end

% find the frequency of interest
fbin = nearest(freq.freq, cfg.frequency);

if isfield(freq, 'powspctrm') && isfield(freq, 'crsspctrm')
  % use the power and cross spectrum and construct a square matrix

  % find the index of each sensor channel into powspctrm
  % keep the channel order of the cfg
  [dum, powspctrmindx] = match_str(cfg.channel, freq.label);
  Nchans = length(powspctrmindx);

  % find the index of each sensor channel combination into crsspctrm
  % keep the channel order of the cfg
  crsspctrmindx = zeros(Nchans);
  for sgncmblop=1:size(freq.labelcmb,1)
    ch1 = find(strcmp(cfg.channel, freq.labelcmb(sgncmblop,1)));
    ch2 = find(strcmp(cfg.channel, freq.labelcmb(sgncmblop,2)));
    if ~isempty(ch1) && ~isempty(ch2)
      % this square matrix contains the indices into the signal combinations
      crsspctrmindx(ch1,ch2) = sgncmblop;
    end
  end

  % this complex rearrangement of channel indices transforms the CSDs into a square matrix
  if strcmp(freq.dimord, 'chan_freq')
    % FIXME this fails in case dimord=rpt_chan_freq and only 1 trial
    Cf = complex(nan*zeros(Nchans,Nchans));
    % first use the complex conjugate for all reversed signal combinations
    Cf(find(crsspctrmindx)) = freq.crsspctrm(crsspctrmindx(find(crsspctrmindx)), fbin);
    Cf = ctranspose(Cf);
    % and then get get the csd for all signal combinations
    Cf(find(crsspctrmindx)) = freq.crsspctrm(crsspctrmindx(find(crsspctrmindx)), fbin);
    % put the power on the diagonal
    Cf(find(eye(Nchans))) = freq.powspctrm(powspctrmindx, fbin);
  else
    Cf  = complex(nan*zeros(Ntrials,Nchans,Nchans));
    tmp = complex(nan*zeros(Nchans,Nchans));
    for trial=1:Ntrials
      % first use the complex conjugate for all signal combinations reversed
      tmp(find(crsspctrmindx)) = freq.crsspctrm(trial, crsspctrmindx(find(crsspctrmindx)), fbin);
      tmp = ctranspose(tmp);
      % and then get get the csd for all signal combinations
      tmp(find(crsspctrmindx)) = freq.crsspctrm(trial, crsspctrmindx(find(crsspctrmindx)), fbin);
      % put the power on the diagonal
      tmp(find(eye(Nchans))) = freq.powspctrm(trial, powspctrmindx, fbin);
      Cf(trial,:,:) = tmp;
    end
  end

  % do a sanity check on the cross-spectral-density matrix
  if any(isnan(Cf(:)))
    error('The cross-spectral-density matrix is not complete');
  end

  if isfield(cfg, 'refchan') && ~isempty(cfg.refchan)
    % contruct the cross-spectral-density vector of the reference channel with all MEG channels
    tmpindx = match_str(freq.labelcmb(:,1), cfg.refchan);
    refindx = match_str(freq.labelcmb(tmpindx,2), cfg.channel);
    refindx = tmpindx(refindx);
    flipref = 0;
    if isempty(refindx)
      % first look in the second column, then in the first
      tmpindx = match_str(freq.labelcmb(:,2), cfg.refchan);
      refindx = match_str(freq.labelcmb(tmpindx,1), cfg.channel);
      refindx = tmpindx(refindx);
      flipref = 1;
    end
    if length(refindx)~=Nchans
      error('The cross-spectral-density with the reference channel is not complete');
    end
    if Ntrials==1
      Cr = freq.crsspctrm(refindx, fbin);
    else
      for trial=1:Ntrials
        Cr(trial,:) = freq.crsspctrm(trial, refindx, fbin);
      end
    end
    if flipref
      Cr = conj(Cr);
    end
    % obtain the power of the reference channel
    refindx = match_str(freq.label, cfg.refchan);
    if length(refindx)<1
      error('The reference channel was not found in powspctrm');
    elseif length(refindx)>1
      error('Multiple occurences of the reference channel found in powspctrm');
    end
    if Ntrials==1
      Pr = freq.powspctrm(refindx, fbin);
    else
      for trial=1:Ntrials
        Pr(trial) = freq.powspctrm(trial, refindx, fbin);
      end
      Pr = Pr(:);	% ensure that the first dimension contains the trials
    end
  end

  if strcmp(cfg.dicsfix, 'yes')
    Cr = conj(Cr);
  end

else
  fprintf('computing cross-spectrum from fourier\n');
  [dum, powspctrmindx] = match_str(cfg.channel, freq.label);
  % use the fourier spectrum to compute the complete cross spectrum matrix
  Nchans = length(powspctrmindx);
  if strcmp(freq.dimord, 'chan_freq')
    error('incompatible dimord for computing CSD matrix from fourier');
  elseif strcmp(freq.dimord, 'rpt_chan_freq')
    error('incompatible dimord for computing CSD matrix from fourier');
  elseif strcmp(freq.dimord, 'rpttap_chan_freq'),
    if quickflag,
      Ntrials = 1;
    end
    Cf = zeros(Ntrials,Nchans,Nchans);
    refindx = match_str(freq.label, cfg.refchan);
    if ~isempty(refindx)
      Cr = zeros(Ntrials,Nchans,1);
      Pr = zeros(Ntrials,1,1);
    end

    if quickflag,
      ntap = sum(freq.cumtapcnt);
      dat  = transpose(freq.fourierspctrm(:, powspctrmindx, fbin));
      Cf(1,:,:) = (dat * ctranspose(dat)) ./ ntap;
      if ~isempty(refindx)
        ref = transpose(freq.fourierspctrm(:, refindx, fbin));
	Cr(1,:,1) = dat * ctranspose(ref) ./ ntap;
	Pr(1,1,1) = ref * ctranspose(ref) ./ ntap;
      end
    else
      freq.cumtapcnt = freq.cumtapcnt(:)';
      for k=1:Ntrials
        tapbeg = 1 + sum([0 freq.cumtapcnt(1:(k-1))]);
        tapend =     sum([0 freq.cumtapcnt(1:(k  ))]);
        ntap = freq.cumtapcnt(k);
        dat  = transpose(freq.fourierspctrm(tapbeg:tapend, powspctrmindx, fbin));
        Cf(k,:,:) = (dat * ctranspose(dat)) ./ ntap;
        if ~isempty(refindx)
          ref = transpose(freq.fourierspctrm(tapbeg:tapend, refindx, fbin));
          Cr(k,:,1) = dat * ctranspose(ref) ./ ntap;
          Pr(k,1,1) = ref * ctranspose(ref) ./ ntap;
        end
      end
    end
  else
    error('unsupported dimord for fourier cross-spectrum computation');
  end
end

% update the configuration, so that the calling function exactly knows what was selected
if ~isempty(tbin),
  % a single latency was selected in the freq structure
  cfg.latency = freq.time;
else
  cfg.latency = [];
end
cfg.frequency = freq.freq(fbin);
cfg.channel   = freq.label(powspctrmindx);

