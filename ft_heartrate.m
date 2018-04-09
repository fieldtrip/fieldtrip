function [dataout] = ft_heartrate(cfg, datain)

% FT_HEARTRATE estimates the heartrate from a continuous PPG or ECG channel. It
% returns a new data structure with a continuous representation of the heartrate in
% beats per minute.
%
% Use as
%   dataout = ft_analyze_heartrate(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING.
%
% The configuration structure has the following options
%   cfg.channel          = selected channel for processing, see FT_CHANNELSELECTION
%   cfg.envelopewindow   = scalar, time in seconds
%   cfg.medianwindow     = integer, number of heartbeats
%   cfg.threshold        = scalar, between 0 and 1
%   cfg.feedback         = 'yes' or 'no'
% The input data can be preprocessed on the fly using
%   cfg.preproc.bpfilter = 'yes' or 'no'
%   cfg.preproc.bpfreq   = [low high], filter frequency in Hz
%
% See also FT_ELECTRODERMALACTIVITY, FT_HEADMOVEMENT, FT_REGRESSCONFOUND

% Copyright (C) 2018, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% check if the input data is valid for this function, the input data must be raw
datain = ft_checkdata(datain, 'datatype', 'raw', 'feedback', 'yes');

% set the default options
cfg.channel          = ft_getopt(cfg, 'channel', {});
cfg.envelopewindow   = ft_getopt(cfg, 'envelopewindow', 10);  % in seconds
cfg.medianwindow     = ft_getopt(cfg, 'medianwindow', 120);   % in heartbeats
cfg.threshold        = ft_getopt(cfg, 'threshold', 0.4);      % between 0 and 1
cfg.feedback         = ft_getopt(cfg, 'feedback', 'yes');
cfg.preproc          = ft_getopt(cfg, 'preproc', []);
cfg.preproc.bpfilter = ft_getopt(cfg.preproc, 'bpfilter', 'yes');
cfg.preproc.bpfreq   = ft_getopt(cfg.preproc, 'bpfreq', [0.2 20]);    % in Hz

% copy some of the fields over to the new data structure
dataout = keepfields(datain, {'time', 'fsample', 'sampleinfo', 'trialinfo'});
dataout.label = {'heartrate'};
dataout.trial = {};  % this is to be determined in the main code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.channel = ft_channelselection(cfg.channel, datain.label);
assert(numel(cfg.channel)==1, 'you should specify exactly one channel');

chansel = strcmp(datain.label, cfg.channel{1});
fsample = datain.fsample;
envelopewindow = round(cfg.envelopewindow*fsample); % in samples

for trllop=1:numel(datain.trial)
  dat   = datain.trial{trllop}(chansel,:);
  label = datain.label(chansel);
  time  = datain.time{trllop};
  
  [yupper,ylower] = envelope(dat, envelopewindow, 'rms');
  
  if istrue(cfg.feedback)
    figure
    subplot(4,1,1)
    hold on
    plot(time, dat)
    plot(time, yupper, 'g');
    plot(time, ylower, 'g');
    xlim([min(time) max(time)])
    xlabel('time (s)');
    title(sprintf('original, trial %d', trllop))
  end
  
  if ~isempty(cfg.preproc)
    % apply the preprocessing to teh selected channel
    [dat, label, time, cfg.preproc] = preproc(dat, label, time, cfg.preproc, 0, 0);
  end
  [yupper,ylower] = envelope(dat, envelopewindow, 'rms');
  
  if istrue(cfg.feedback)
    subplot(4,1,2)
    hold on
    plot(time, dat)
    plot(time, yupper, 'g');
    plot(time, ylower, 'g');
    xlim([min(time) max(time)])
    xlabel('time (s)');
    title('filtered')
  end
  
  dat = (dat - ylower) ./ (yupper - ylower);
  [yupper,ylower] = envelope(dat, envelopewindow, 'rms');
  
  if istrue(cfg.feedback)
    subplot(4,1,3)
    hold on
    plot(time, dat)
    plot(time, yupper, 'g');
    plot(time, ylower, 'g');
    xlim([min(time) max(time)])
    xlabel('time (s)');
    title('locally rescaled')
  end
  
  % find the sample numbers where the filtered value increases above the threshold
  sample = find(diff([0 dat>cfg.threshold])>0);
  ibi = diff(sample)/fsample;
  ibi = ft_preproc_medianfilter(ibi, cfg.medianwindow);
  bpm = 60./ibi;
  
  % compute the sample number right in between beat N and N+1
  sample = round((sample(1:end-1) + sample(2:end))/2);
  
  if istrue(cfg.feedback)
    subplot(4,1,4)
    plot(time(sample), bpm)
    ylim([0 160])
    xlim([min(time) max(time)])
    xlabel('time (s)');
    ylabel('rate (bpm)');
  end
  
  % interpolate it onto the time points of the original data set
  dataout.trial{trllop} = interp1(sample./fsample, bpm, time);
end % for trllop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout
