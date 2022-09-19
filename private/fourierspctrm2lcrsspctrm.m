function [freqout] = fourierspctrm2lcrsspctrm(freq, varargin)

% FOURIERSPCTRM2LCRSSPCTRM is a helper function that converts the input fourierspctrm
% into a lagged crsspctrm, to enable computation of lagged coherence as described in
% the publication referenced below. It is used in ft_connectivityanalysis in order to
% reorganize the input data.
%
% The input data should be organised in a structure as obtained from the
% FT_FREQANALYSIS function (freq), such that freq contains the fields 'fourierspctrm'
% and 'time'. The timepoints must be chosen such that the desired cfg.lag/cfg.foi
% (lag in seconds) is an integer multiple of the time resolution in freq.
%
% Options come in key-value pairs, and may contain
%   lag = scalar (or vector) of time shifts, expressed in units of time
%      We recommend users to choose cfg.lag such that it is larger or equal
%      to the width of the wavelet used for each Fourier transform in ft_freqanalysis
%   timeresolved = 'yes' or 'no' (default='no'). If set to yes, lagged
%      coherence is calculated separately for each pair of timepoints that
%      is separated by lag
%   channelcmb   =  Mx2 cell-array with selection of channel pairs,
%      see ft_channelcombination, default = {'all' 'all'};
%
% When this measure is used for your publication, please cite:
%   Fransen, Anne M. M, Van Ede, Freek, Maris, Eric (2015) Identifying
%   oscillations on the basis of rhythmicity. NeuroImage 118: 256-267.
% You may also want to acknowledge the fact that J.M. Schoffelen has
% written the actual implementation.

% Copyright (C) 2019, Jan-Mathijs Schoffelen, DCCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
% FieldTrip is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% FieldTrip is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%

lags         = ft_getopt(varargin,        'lags',         []);
timeresolved = istrue(ft_getopt(varargin, 'timeresolved', 'no'));
channelcmb   = ft_getopt(varargin,        'channelcmb',   {'all' 'all'});

% check if the input data is valid for this function
if ~isfield(freq, 'time')
  ft_error('this function requires frequency data with a time axis')
end

% deal with the lag, if it's empty, throw an error if number of lags>1 and timeresolved is true
if isempty(lags)
  lags = mean(diff(freq.time));
end
if numel(lags)>1 && timeresolved
  error('in case of a timeresolved estimate, only a single time lag is allowed');
end

% deal with the channelcmb -> not needed: auto combinations are already
% handled by the invoking function
%channelcmb = ft_channelcombination(channelcmb, freq.label, 1);

% determine the corresponding indices of all channel combinations
[dummy, chancmbind(:,1)] = match_str(channelcmb(:,1), freq.label);
[dummy, chancmbind(:,2)] = match_str(channelcmb(:,2), freq.label);
clear dummy

ntime = numel(freq.time);

% convert the lags, expressed in units of time, to samples
lagsorig  = lags;
sellags   = lags<=freq.time(end)-freq.time(1);
lags_indx = round(lags(sellags)./mean(diff(freq.time)));
fprintf('\nrequested lag(s): %0.3d\t', lags);
lags = lags_indx.*mean(diff(freq.time));
fprintf('\nused lags(s): %0.3d\t', lags);

label = freq.label;
label_lagged = freq.label;
for k = 1:numel(label)
  label_lagged{k,1} = sprintf('%s_lagged',label{k});
end
tmpchannelcmb = [label(chancmbind(:,1)) label_lagged(chancmbind(:,2));label label;label_lagged label_lagged];

for k = 1:numel(lags_indx)
  indx1   = 1:(ntime-lags_indx(k));
  indx2   = (1+lags_indx(k)):ntime;

  tmpfreq = freq;
  tmpfreq.fourierspctrm = cat(2, tmpfreq.fourierspctrm(:,:,:,indx1), ...
                                 tmpfreq.fourierspctrm(:,:,:,indx2));
  tmpfreq.time  = freq.time(indx1);
  tmpfreq.label = cat(1, label, label_lagged);
  tmpfreq       = ft_checkdata(tmpfreq, 'cmbstyle', 'sparse', 'channelcmb', tmpchannelcmb);

  if timeresolved
    error('to do');
  else
    % convert the data matrix into a 'rpt_chancmb_freq', accounting for
    % time points with nans may be needed
    ntimex = numel(tmpfreq.time);
    if ~contains(tmpfreq.dimord, 'rpt')
      nrpt = 1;
      permutevec = [3 1 2];
    else
      nrpt  = size(tmpfreq.crsspctrm,1);
      permutevec = [4 1 2 3];
    end
    ncmb  = size(tmpfreq.labelcmb,1);
    nfreq = numel(tmpfreq.freq);
    tmpfreq.lcrsspctrm = reshape(permute(tmpfreq.crsspctrm, permutevec),[nrpt*ntimex ncmb nfreq]);
    tmpfreq = removefields(tmpfreq, {'trialinfo'});
    tmpfreq.dimord = 'rpt_chancmb_freq_time';

    if k==1
      tmpfreq.cumtapcnt  = repmat(tmpfreq.cumtapcnt, [ntimex 1]);

      freqout = removefields(tmpfreq, {'time' 'crsspctrm'});
      freqout.time = zeros(1,0);
      freqout.lcrsspctrm(:,:,:,2:numel(lags_indx)) = nan;
    else
      ntmp = size(tmpfreq.lcrsspctrm,1);
      freqout.lcrsspctrm(1:ntmp,:,:,k) = tmpfreq.lcrsspctrm;
    end
    freqout.time      = cat(2, freqout.time, lags(k)); % use the 'time' as fieldname to avoid bookkeeping crashes
  end
end

% use only those slices that contain data for all channel pairs, to ensure
% the same degrees of freedom for denominator and numerator: note that the
% degrees of freedom can go down as a function of time lag

freqout.lcrsspctrm(repmat(~all(isfinite(freqout.lcrsspctrm),2), [1 ncmb 1 1]))=nan;

% check for a discrepancy between the requested lags and the output data
if ~all(sellags)
  siz = size(freqout.lcrsspctrm);
  if numel(freqout.freq)==1
    siz(3) = 1;
  end
  siz(4) = numel(lagsorig);
  tmp = nan(siz);
  tmp(:,:,:,sellags) = freqout.lcrsspctrm;
  freqout.lcrsspctrm = tmp;
  clear tmp;
  time = lagsorig;
  time(sellags) = freqout.time;
  freqout.time  = time;
end
