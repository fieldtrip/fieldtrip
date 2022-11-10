function [data] = ft_dipolesimulation(cfg)

% FT_DIPOLESIMULATION simulates channel-level time-series data that consists of the
% the spatial distribution of the the field or potential of one or multiple dipoles.
%
% Use as
%   data = ft_dipolesimulation(cfg)
% which will return a raw data structure that resembles the output of
% FT_PREPROCESSING.
%
% The dipoles position and orientation have to be specified with
%   cfg.sourcemodel.pos        = [Rx Ry Rz] (size Nx3)
%   cfg.sourcemodel.mom        = [Qx Qy Qz] (size 3xN)
%   cfg.sourcemodel.unit       = string, can be 'mm', 'cm', 'm' (default is automatic)
%
% The timecourse of the dipole activity is given as a cell-array with one
% dipole signal per trial
%   cfg.sourcemodel.signal     = cell-array with one dipole signal per trial
% or by specifying the parameters of a sine-wave signal
%   cfg.sourcemodel.frequency  =   in Hz
%   cfg.sourcemodel.phase      =   in radians
%   cfg.sourcemodel.amplitude  =   per dipole
%
% The number of trials and the time axes of the trials can be specified by
%   cfg.fsample     = simulated sample frequency (default = 1000)
%   cfg.trllen      = length of simulated trials in seconds (default = 1)
%   cfg.numtrl      = number of simulated trials (default = 10)
%   cfg.baseline    = number (default = 0.3)
% or by
%   cfg.time        = cell-array with one time axis per trial, for example obtained from an existing dataset
%
% Random white noise can be added to the data in each trial, either by
% specifying an absolute or a relative noise level
%   cfg.relnoise    = add noise with level relative to data signal
%   cfg.absnoise    = add noise with absolute level
%   cfg.randomseed  = 'yes' or a number or vector with the seed value (default = 'yes')
%
% Optional input arguments are
%   cfg.channel    = Nx1 cell-array with selection of channels (default = 'all'),
%                    see FT_CHANNELSELECTION for details
%   cfg.dipoleunit = units for dipole amplitude (default nA*m)
%   cfg.chanunit   = Nx1 cell-array with units for the channel data
%
% Optionally, you can modify the leadfields by reducing the rank, i.e. remove the weakest orientation
%   cfg.reducerank    = 'no', or number (default = 3 for EEG, 2 for MEG)
%   cfg.backproject   = 'yes' or 'no',  determines when reducerank is applied whether the
%                       lower rank leadfield is projected back onto the original linear
%                       subspace, or not (default = 'yes')
%
% The volume conduction model of the head should be specified as
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%
% The EEG or MEG sensor positions should be specified as
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition or filename, see FT_READ_SENS
%
% See also FT_SOURCEANALYSIS, FT_DIPOLEFITTING, FT_TIMELOCKSIMULATION,
% FT_FREQSIMULATION, FT_CONNECTIVITYSIMULATION

% Undocumented local options
% cfg.feedback
% cfg.previous
% cfg.version

% Copyright (C) 2004, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance
ft_preamble randomseed

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamed',    {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed',    {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed',    {'optofile', 'opto'});
cfg = ft_checkconfig(cfg, 'renamed',    {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'dip',    'sourcemodel'});

% for consistency with FT_TIMELOCKSIMULUATION and FT_FREQSIMULATION
cfg = ft_checkconfig(cfg, 'createsubcfg', 'sourcemodel');
cfg = ft_checkconfig(cfg, 'renamed', {'ntrials', 'numtrl'});
cfg = ft_checkconfig(cfg, 'renamed', {'triallength', 'trllen'});

% set the defaults
cfg.sourcemodel         = ft_getopt(cfg, 'sourcemodel', []);
cfg.sourcemodel.pos     = ft_getopt(cfg.sourcemodel, 'pos', [-5 0 15]);
cfg.sourcemodel.mom     = ft_getopt(cfg.sourcemodel, 'mom', [1 0 0]');
cfg.sourcemodel.time    = ft_getopt(cfg.sourcemodel, 'time', {});
cfg.sourcemodel.signal  = ft_getopt(cfg.sourcemodel, 'signal', {});
cfg.fsample             = ft_getopt(cfg, 'fsample', 250);
cfg.relnoise            = ft_getopt(cfg, 'relnoise', 0);
cfg.absnoise            = ft_getopt(cfg, 'absnoise', 0);
cfg.feedback            = ft_getopt(cfg, 'feedback', 'text');
cfg.channel             = ft_getopt(cfg, 'channel',  'all');
cfg.dipoleunit          = ft_getopt(cfg, 'dipoleunit', 'nA*m');
cfg.chanunit            = ft_getopt(cfg, 'chanunit', {});

% collect and preprocess the electrodes/gradiometer and head model
% this will also update cfg.channel to match the electrodes/gradiometers
[headmodel, sens, cfg] = prepare_headmodel(cfg, []);

% construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD
leadfieldopt = {};
leadfieldopt = ft_setopt(leadfieldopt, 'reducerank',     ft_getopt(cfg, 'reducerank'));
leadfieldopt = ft_setopt(leadfieldopt, 'backproject',    ft_getopt(cfg, 'backproject'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(cfg, 'normalize'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(cfg, 'normalizeparam'));
leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(cfg, 'weight'));

cfg.sourcemodel = fixdipole(cfg.sourcemodel);
Ndipoles = size(cfg.sourcemodel.pos,1);

% in case no time or signal was given, set some additional defaults
if ~isempty(cfg.sourcemodel.time) && ~isempty(cfg.sourcemodel.signal)
  assert(length(cfg.sourcemodel.signal)==length(cfg.sourcemodel.time)); % these must match
  cfg.numtrl    = length(cfg.sourcemodel.time);
  cfg.fsample   = 1/mean(diff(cfg.sourcemodel.time{1}));  % determine from time-axis
  cfg.trllen    = length(cfg.sourcemodel.time{1})/cfg.fsample;
  cfg.baseline  = -cfg.sourcemodel.time{1}(1);
elseif ~isempty(cfg.sourcemodel.time)
  cfg.numtrl    = length(cfg.sourcemodel.time);
  cfg.fsample   = 1/mean(diff(cfg.sourcemodel.time{1}));  % determine from time-axis
  cfg.trllen    = length(cfg.sourcemodel.time{1})/cfg.fsample;
  cfg.baseline  = -cfg.sourcemodel.time{1}(1);
elseif ~isempty(cfg.sourcemodel.signal)
  cfg.numtrl    = length(cfg.sourcemodel.signal);
  cfg.fsample   = ft_getopt(cfg, 'fsample', 1000);
  cfg.trllen    = length(cfg.sourcemodel.signal{1})/cfg.fsample;
  cfg.baseline  = ft_getopt(cfg, 'baseline', 0);
else
  cfg.numtrl    = ft_getopt(cfg, 'numtrl', 10);
  cfg.fsample   = ft_getopt(cfg, 'fsample', 1000);
  cfg.trllen    = ft_getopt(cfg, 'trllen', 1);
  cfg.baseline  = ft_getopt(cfg, 'baseline', 0);
end

% no signal was given, set some additional defaults
if isempty(cfg.sourcemodel.signal)
  cfg.sourcemodel.frequency   = ft_getopt(cfg, 'frequency', ones(Ndipoles,1)*10);
  cfg.sourcemodel.phase       = ft_getopt(cfg, 'phase', zeros(Ndipoles,1));
  cfg.sourcemodel.amplitude   = ft_getopt(cfg, 'amplitude', ones(Ndipoles,1));
end

if isfield(cfg.sourcemodel, 'frequency')
  % this should be a column vector
  cfg.sourcemodel.frequency = cfg.sourcemodel.frequency(:);
end

if isfield(cfg.sourcemodel, 'phase')
  % this should be a column vector
  cfg.sourcemodel.phase = cfg.sourcemodel.phase(:);
end

if ~isempty(cfg.sourcemodel.time)
  % use the user-supplied time vectors
  diptime = cfg.sourcemodel.time;
else
  % construct a time axis for every trial
  nsample = round(cfg.trllen*cfg.fsample);
  diptime = cell(1, cfg.numtrl);
  for iTr = 1:cfg.numtrl
    diptime{iTr} = (((1:nsample)-1)/cfg.fsample) - cfg.baseline;
  end
end

if ~isempty(cfg.sourcemodel.signal)
  % use the user-supplied signal for the dipoles
  dipsignal = cfg.sourcemodel.signal;
else
  dipsignal = cell(1, cfg.numtrl);
  for iTr = 1:cfg.numtrl
    % compute a cosine signal with the desired frequency, phase and amplitude for each dipole
    for i=1:Ndipoles
      dipsignal{iTr}(i,:) = cos(cfg.sourcemodel.frequency(i)*diptime{iTr}*2*pi + cfg.sourcemodel.phase(i)) * cfg.sourcemodel.amplitude(i);
    end
  end
end

dippos = cfg.sourcemodel.pos;
dipmom = cfg.sourcemodel.mom;

if ~iscell(dipmom)
  dipmom = {dipmom};
end

if ~iscell(dippos)
  dippos = {dippos};
end

if length(dippos)==1
  dippos = repmat(dippos, 1, cfg.numtrl);
elseif length(dippos)~=cfg.numtrl
  ft_error('incorrect number of trials specified in the dipole position');
end

if length(dipmom)==1
  dipmom = repmat(dipmom, 1, cfg.numtrl);
elseif length(dipmom)~=cfg.numtrl
  ft_error('incorrect number of trials specified in the dipole moment');
end

data.time  = diptime;
data.trial = {};

ft_progress('init', cfg.feedback, 'computing data data');
for trial=1:cfg.numtrl
  ft_progress(trial/cfg.numtrl, 'computing data data for trial %d\n', trial);
  if numel(cfg.chanunit) == numel(cfg.channel)
    lf = ft_compute_leadfield(dippos{trial}, sens, headmodel, 'dipoleunit', cfg.dipoleunit, 'chanunit', cfg.chanunit, leadfieldopt{:});
  else
    lf = ft_compute_leadfield(dippos{trial}, sens, headmodel, leadfieldopt{:});
  end
  nsamples = size(dipsignal{trial},2);
  nchannels = size(lf,1);
  data.trial{trial} = zeros(nchannels,nsamples);
  for i = 1:3
    data.trial{trial} = data.trial{trial} + ...
      lf(:,i:3:end) * (repmat(dipmom{trial}(i:3:end),1,nsamples) .* dipsignal{trial});
  end
end
ft_progress('close');

if ft_senstype(sens, 'meg')
  data.grad = sens;
elseif ft_senstype(sens, 'eeg')
  data.elec = sens;
end

% determine RMS value of data data
ss = 0;
sc = 0;
for trial=1:cfg.numtrl
  ss = ss + sum(data.trial{trial}(:).^2);
  sc = sc + length(data.trial{trial}(:));
end
rms = sqrt(ss/sc);
ft_info('RMS value of data data is %g\n', rms);

% add noise to the data data
for trial=1:cfg.numtrl
  relnoise = randn(size(data.trial{trial})) * cfg.relnoise * rms;
  absnoise = randn(size(data.trial{trial})) * cfg.absnoise;
  data.trial{trial} = data.trial{trial} + relnoise + absnoise;
end

data.fsample = cfg.fsample;
data.label   = sens.label;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble randomseed
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
