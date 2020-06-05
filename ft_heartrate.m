function [dataout] = ft_heartrate(cfg, datain)

% FT_HEARTRATE estimates the heart rate from a continuous PPG or ECG channel. It
% returns a new data structure with a continuous representation of the heartrate in
% beats per minute.
%
% Use as
%   dataout = ft_heartrate(cfg, data)
% where the input data is a structure as obtained from FT_PREPROCESSING.
%
% The configuration structure has the following general options
%   cfg.channel          = selected channel for processing, see FT_CHANNELSELECTION
%   cfg.feedback         = 'yes' or 'no'
%   cfg.method           = string representing the method for heart beat detection
%                          'findpeaks'  filtering and normalization, followed by FINDPEAKS (default)
%                          'pantompkin' implementation of the Pan-Tompkin algorithm for ECG beat detection
%
% For the 'findpeaks' method the following additional options can be specified
%   cfg.envelopewindow   = scalar, time in seconds (default = 10)
%   cfg.peakseparation   = scalar, time in seconds
%   cfg.threshold        = scalar, usually between 0 and 1 (default = 0.4)
%   cfg.flipsignal       = 'yes' or 'no', whether to flip the polarity of the signal (default is automatic)
% and the data can be preprocessed on the fly using
%   cfg.preproc.bpfilter = 'yes' or 'no'
%   cfg.preproc.bpfreq   = [low high], filter frequency in Hz
% This implementation performs some filtering and amplitude normalization, followed
% by the FINDPEAKS function. It works both for ECG as for PPG signals.
%
% For the 'pantompkin` method there are no additional options. This implements
% - J Pan, W J Tompkins, "A Real-Time QRS Detection Algorithm", IEEE Trans Biomed Eng, 1985. https://doi.org/10.1109/tbme.1985.325532
% - H Sedghamiz, "Matlab Implementation of Pan Tompkins ECG QRS detector". https://doi.org/10.13140/RG.2.2.14202.59841
%
% See also FT_ELECTRODERMALACTIVITY, FT_HEADMOVEMENT, FT_REGRESSCONFOUND

% Copyright (C) 2018-2020, Robert Oostenveld & Helena Cockx
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

% ensure that users with old scripts are aware of changes
cfg = ft_checkconfig(cfg, 'forbidden', 'medianwindow');

% set the default options
cfg.channel          = ft_getopt(cfg, 'channel', 'all');
cfg.method           = ft_getopt(cfg, 'method', 'findpeaks');
cfg.envelopewindow   = ft_getopt(cfg, 'envelopewindow', 10);  % in seconds
cfg.peakseparation   = ft_getopt(cfg, 'peakseparation', []);  % in seconds
cfg.threshold        = ft_getopt(cfg, 'threshold', 0.4);      % between 0 and 1
cfg.feedback         = ft_getopt(cfg, 'feedback', 'yes');
cfg.preproc          = ft_getopt(cfg, 'preproc', []);
cfg.flipsignal       = ft_getopt(cfg, 'flipsignal', []);

% the expected rate is around 80 bpm, which means 80/60=1.33 Hz
cfg.preproc.bpfilter    = ft_getopt(cfg.preproc, 'bpfilter', 'yes');
cfg.preproc.bpfilttype  = ft_getopt(cfg.preproc, 'bpfilttype', 'but');
cfg.preproc.bpfiltdir   = ft_getopt(cfg.preproc, 'bpfiltdir', 'twopass');
cfg.preproc.bpfiltord   = ft_getopt(cfg.preproc, 'bpfiltord', 2);
cfg.preproc.bpfreq      = ft_getopt(cfg.preproc, 'bpfreq', [1/3 10] * 1.33);  % in Hz

% copy some of the fields over to the new data structure
dataout = keepfields(datain, {'time', 'fsample', 'sampleinfo', 'trialinfo'});
dataout.label = {'heartrate', 'heartbeatphase', 'heartbeatonset'};
dataout.trial = {};  % this is to be determined in the main code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.channel = ft_channelselection(cfg.channel, datain.label);
assert(numel(cfg.channel)==1, 'you should specify exactly one channel');

chansel = strcmp(datain.label, cfg.channel{1});
fsample = datain.fsample;

switch cfg.method
  case 'findpeaks'
    for trllop=1:numel(datain.trial)
      dat   = datain.trial{trllop}(chansel,:);
      label = datain.label(chansel);
      time  = datain.time{trllop};
      
      if isempty(cfg.flipsignal)
        if skewness(dat)<0
          cfg.flipsignal = 'yes';
        else
          cfg.flipsignal = 'no';
        end
      end
      
      if istrue(cfg.flipsignal)
        ft_notice('flipping signal polarity');
        dat = -dat;
      end
      
      if ~isempty(cfg.peakseparation)
        [yupper,ylower] = envelope(dat, round(cfg.peakseparation*fsample), 'peaks');
      elseif ~isempty(cfg.envelopewindow)
        [yupper,ylower] = envelope(dat, round(cfg.envelopewindow*fsample), 'rms');
      end
      
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
        % apply the preprocessing to the selected channel
        [dat, label, time, cfg.preproc] = preproc(dat, label, time, cfg.preproc, 0, 0);
      end
      
      if ~isempty(cfg.peakseparation)
        [yupper,ylower] = envelope(dat, round(cfg.peakseparation*fsample), 'peaks');
      elseif ~isempty(cfg.envelopewindow)
        [yupper,ylower] = envelope(dat, round(cfg.envelopewindow*fsample), 'rms');
      end
      
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
      
      if ~isempty(cfg.peakseparation)
        [yupper,ylower] = envelope(dat, round(cfg.peakseparation*fsample), 'peaks');
      elseif ~isempty(cfg.envelopewindow)
        [yupper,ylower] = envelope(dat, round(cfg.envelopewindow*fsample), 'rms');
      end
      
      % find the sample numbers where the filtered value increases above the threshold
      [vals, peaks] = findpeaks(dat, 'MinPeakHeight', cfg.threshold);
      
      if istrue(cfg.feedback)
        subplot(4,1,3)
        hold on
        plot(time, dat)
        plot(time, yupper, 'g');
        plot(time, ylower, 'g');
        plot(time(peaks), vals, 'r*');
        xlim([min(time) max(time)])
        xlabel('time (s)');
        title('locally rescaled')
      end
      
      % construct a continuous channel with the rate and the phase
      [rate, phase, tmp] = discr2ctu(peaks, size(dat), fsample);
      
      % add the continuous channels to the output structure
      dataout.trial{trllop} = [rate; phase; tmp];
      
      if istrue(cfg.feedback)
        subplot(4,1,4)
        plot(time, rate)
        ylim([0 160])
        xlim([min(time) max(time)])
        xlabel('time (s)');
        ylabel('rate (bpm)');
      end
      
      ft_info('heart rate in trial %d: mean=%.1f, min=%.1f, max=%.1f\n', trllop, nanmean(rate), nanmin(rate), nanmax(rate));
      
    end % for trllop
    
  case 'pantompkin'
    ft_hastoolbox('fileexchange', 1);
    for trllop=1:numel(datain.trial)
      dat   = datain.trial{trllop}(chansel,:);
      label = datain.label(chansel);
      time  = datain.time{trllop};
      
      % pan-tompkin algorithm
      [vals, peaks, delay] = pan_tompkin(dat, fsample, istrue(cfg.feedback));
      
      % construct a continuous channel with the rate and the phase
      [rate, phase, tmp] = discr2ctu(peaks, size(dat), fsample);
      
      % add the continuous channels to the output structure
      dataout.trial{trllop} = [rate; phase; tmp];
    end
    
  otherwise
    ft_error('unsupported method %s', cfg.method);
    
end % switch method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rate, phase, tmp] = discr2ctu(peaks, n, fsample)
rate  = nan(n);
phase = nan(n);
for i=1:length(peaks)-1
  begsample = peaks(i);
  endsample = peaks(i+1);
  rate(begsample:endsample)  = 60 * fsample/(endsample-begsample); % in bpm
  phase(begsample:endsample) = linspace(-pi, pi, (endsample-begsample+1));
end
% also construct a boolean channel with a pulse at the beat onset
tmp = zeros(n);
tmp(peaks) = 1;
