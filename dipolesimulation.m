function [simulated] = dipolesimulation(cfg)

% DIPOLESIMULATION computes the field or potential of a simulated dipole
% and returns a datastructure identical to the PREPROCESSING function.
%
% Use as
%   data = dipolesimulation(cfg)
%
% You should specify the volume conductor model with
%   cfg.hdmfile       = string, file containing the volume conduction model
% or alternatively
%   cfg.vol           = structure with volume conduction model
% If the sensor information is not contained in the data itself you should 
% also specify the sensor information using
%   cfg.gradfile      = string, file containing the gradiometer definition
%   cfg.elecfile      = string, file containing the electrode definition
% or alternatively
%   cfg.grad          = structure with gradiometer definition
%   cfg.elec          = structure with electrode definition
%
% optionally
%   cfg.channel    = Nx1 cell-array with selection of channels (default = 'all'),
%                    see CHANNELSELECTION for details
%
% The dipoles position and orientation have to be specified with
%   cfg.dip.pos     = [Rx Ry Rz]
%   cfg.dip.mom     = [Qx Qy Qz]
%
% The timecourse of the dipole activity is given as a single vector or as a
% cell-array with one vectors per trial
%   cfg.dip.signal
% or by specifying a sine-wave signal
%   cfg.dip.frequency    in Hz
%   cfg.dip.phase        in radians
%   cfg.ntrials          number of trials
%   cfg.triallength      time in seconds
%   cfg.fsample          sampling frequency in Hz
%
% Random white noise can be added to the data in each trial, either by 
% specifying an absolute or a relative noise level
%   cfg.relnoise    = add noise with level relative to simulated signal
%   cfg.absnoise    = add noise with absolute level

% Undocumented local options
% cfg.feedback
% cfg.previous
% cfg.version
%
% This function depends on PREPARE_VOL_SENS which has the following options:
% cfg.channel, documented
% cfg.elec, documented
% cfg.elecfile, documented
% cfg.grad, documented
% cfg.gradfile, documented
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: dipolesimulation.m,v $
% Revision 1.24  2009/07/02 15:37:32  roboos
% use fixdipole for consistent dipole structure representation
%
% Revision 1.23  2009/03/23 21:19:51  roboos
% allow different amplitudes for different dipoles
%
% Revision 1.22  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.21  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.20  2008/03/18 12:27:22  roboos
% use senstype helper function to determine whether output should contain grad or elec
%
% Revision 1.19  2007/03/21 12:43:41  roboos
% replaced two "try" statements with if(isfield(..))
%
% Revision 1.18  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.17  2006/07/04 16:05:33  roboos
% removed offset from output data
%
% Revision 1.16  2006/06/20 16:25:58  ingnie
% updated documentation
%
% Revision 1.15  2006/06/13 14:51:38  ingnie
% updated documentation to increase consistency in help of cfg options, added
% default cfg.channel = 'all'
%
% Revision 1.14  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.13  2006/04/10 16:33:46  ingnie
% updated documentation
%
% Revision 1.12  2005/11/10 10:02:14  roboos
% minor change in whitespace and help
%

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'dip'),        cfg.dip = [];             end
if ~isfield(cfg.dip, 'pos'),    cfg.dip.pos = [-5 0 15];  end
if ~isfield(cfg.dip, 'mom'),    cfg.dip.mom = [1 0 0]';   end
if ~isfield(cfg, 'fsample'),    cfg.fsample = 250;        end
if ~isfield(cfg, 'relnoise'),   cfg.relnoise = 0;         end
if ~isfield(cfg, 'absnoise'),   cfg.absnoise = 0;         end
if ~isfield(cfg, 'feedback'),   cfg.feedback = 'text';    end
if ~isfield(cfg, 'channel'),    cfg.channel = 'all';      end

cfg.dip = fixdipole(cfg.dip);
Ndipoles = size(cfg.dip.pos,1);

% prepare the volume conductor and the sensor array
[vol, sens, cfg] = prepare_headmodel(cfg, []);

if ~isfield(cfg, 'ntrials') 
  if isfield(cfg.dip, 'signal')
    cfg.ntrials = length(cfg.dip.signal);
  else
    cfg.ntrials = 20;
  end
end
Ntrials  = cfg.ntrials;

if isfield(cfg.dip, 'frequency')
  % this should be a column vector
  cfg.dip.frequency = cfg.dip.frequency(:);
end

if isfield(cfg.dip, 'phase')
  % this should be a column vector
  cfg.dip.phase = cfg.dip.phase(:);
end

% no signal was given, compute a cosine-wave signal as timcourse for the dipole
if ~isfield(cfg.dip, 'signal')
  % set some additional defaults if neccessary
  if ~isfield(cfg.dip, 'frequency')
    cfg.dip.frequency = ones(Ndipoles,1)*10;
  end
  if ~isfield(cfg.dip, 'phase')
    cfg.dip.phase = zeros(Ndipoles,1);
  end
  if ~isfield(cfg.dip, 'amplitude')
    cfg.dip.amplitude = ones(Ndipoles,1);
  end
  if ~isfield(cfg, 'triallength')
    cfg.triallength = 1;
  end
  % compute a cosine-wave signal wit the desired frequency, phase and amplitude for each dipole
  nsamples = round(cfg.triallength*cfg.fsample);
  time     = (0:(nsamples-1))/cfg.fsample;
  for i=1:Ndipoles
    cfg.dip.signal(i,:) = cos(cfg.dip.frequency(i)*time*2*pi + cfg.dip.phase(i)) * cfg.dip.amplitude(i);
  end
end

% construct the timecourse of the dipole activity for each individual trial
if ~iscell(cfg.dip.signal)
  dipsignal = {};
  time      = {};
  nsamples  = length(cfg.dip.signal);
  for trial=1:Ntrials
    % each trial has the same dipole signal 
    dipsignal{trial} = cfg.dip.signal;
    time{trial} = (0:(nsamples-1))/cfg.fsample;
  end
else
  dipsignal = {};
  time      = {};
  for trial=1:Ntrials
    % each trial has a different dipole signal 
    dipsignal{trial} = cfg.dip.signal{trial};
    time{trial} = (0:(length(dipsignal{trial})-1))/cfg.fsample;
  end
end

dippos    = cfg.dip.pos;
dipmom    = cfg.dip.mom;

if ~iscell(dipmom)
  dipmom = {dipmom};
end

if ~iscell(dippos)
  dippos = {dippos};
end

if length(dippos)==1
  dippos = repmat(dippos, 1, Ntrials);
elseif length(dippos)~=Ntrials
  error('incorrect number of trials specified in the dipole position');
end

if length(dipmom)==1
  dipmom = repmat(dipmom, 1, Ntrials);
elseif length(dipmom)~=Ntrials
  error('incorrect number of trials specified in the dipole moment');
end

simulated.trial  = {};
simulated.time   = {};
progress('init', cfg.feedback, 'computing simulated data');
for trial=1:Ntrials
  progress(trial/Ntrials, 'computing simulated data for trial %d\n', trial);
  lf = compute_leadfield(dippos{trial}, sens, vol);
  simulated.trial{trial}  = lf * dipmom{trial} * dipsignal{trial};
  simulated.time{trial}   = time{trial};
end
progress('close');

if senstype(sens, 'meg')
  simulated.grad = sens;
elseif senstype(sens, 'meg')
  simulated.elec = sens;
end

% determine RMS value of simulated data
ss = 0;
sc = 0;
for trial=1:Ntrials
  ss = ss + sum(simulated.trial{trial}(:).^2);
  sc = sc + length(simulated.trial{trial}(:));
end
rms = sqrt(ss/sc);
fprintf('RMS value of simulated data is %g\n', rms);

% add noise to the simulated data
for trial=1:Ntrials
  relnoise = randn(size(simulated.trial{trial})) * cfg.relnoise * rms;
  absnoise = randn(size(simulated.trial{trial})) * cfg.absnoise;
  simulated.trial{trial} = simulated.trial{trial} + relnoise + absnoise;
end

simulated.fsample = cfg.fsample;
simulated.label   = sens.label;

% add version details to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: dipolesimulation.m,v 1.24 2009/07/02 15:37:32 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
simulated.cfg = cfg;

