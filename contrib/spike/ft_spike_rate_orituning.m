function [Stat] = ft_spike_rate_orituning(cfg, varargin)

% FT_SPIKE_RATE_ORITUNING computes a model of the firing rate as a function
% of orientation or direction.
%
% Use as
%   [stat] = ft_spike_rate_tuning(cfg, rate1, rate2, ... rateN)
%
% The inputs RATE should be the output from FT_SPIKE_RATE. 
%
% Configurations:
%   cfg.stimuli  = should be an 1 x nConditions array of orientations or
%                  directions in radians
%                  varargin{i} corresponds to cfg.stimuli(i)
%   cfg.method   = model to apply, implemented are 'orientation' and 'direction'
%
% Outputs:
%   stat.ang       = mean angle of orientation / direction (1 x nUnits)
%   stat.osi       = orientation selectivity index (Womelsdorf et al., 2012,
%                    PNAS), that is resultant length.
%                    if cfg.method = 'orientation', then orientations are
%                    first projected on the unit circle.
%   stat.di        = direction index, 1 - min/max response

% FIXME: models for contrast etc.

% Copyright (C) 2010, Martin Vinck
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
ft_preamble provenance varargin


% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'stimuli', 'method'});
cfg = ft_checkopt(cfg, 'stimuli', 'doublevector');
cfg = ft_checkopt(cfg, 'method', 'char', {'orientation', 'direction'});

if length(varargin)<2, error('can only compute orituning if multiple inputs are specified'); end 
  
% check whether trials were kept in the rate function
for k = 1:length(varargin)
  try
    varargin{k} = rmfield(varargin{k}, 'trial');
  end
end

Tune        = ft_appendtimelock([],varargin{:});
Tune.avg    = Tune.trial;

% get the unique stimuli and the number of directions, these should match
stimuli     = cfg.stimuli(:);
nConditions = length(varargin);
nStimuli    = length(stimuli);
if nConditions~=nStimuli, error('Length of cfg.stimuli should match number of data inputs'); end

nUnits      = size(Tune.avg,2);

% change the directions so it is a circular variable with range 2*pi
if strcmp(cfg.method,'orientation')
  if (max(stimuli)-min(stimuli))>pi
    error('If cfg.tuningtype is "orientation", cfg.stimuli should have range of pi');
  end
  stimuli = stimuli*2; % convert to make it circular (see Womelsdorf et al. 2012, PNAS).
  
  if (max(stimuli)-min(stimuli))<0.5*pi
    warning('Orientations have a range < 1/2 pi. Are you sure this is correct?. Stats will be biased');
  end
  
  % compute the directionality index
  Stat.di = 1 - min(Tune.avg)./max(Tune.avg);
  
  % transform the data into complex numbers to compute resultant length
  z = exp(1i*stimuli(:)*ones(1,nUnits));
  sumAvg = sum(Tune.avg);
  z = Tune.avg.*z./(sumAvg(ones(nStimuli,1),:));
  Stat.osi  = abs(sum(z));
  
  % make preferred angle modulo 2pi, convert it back to range pi and convert to rad or deg
  prefAngle           = angle(nansum(z));
  prefAngle           = mod(prefAngle,2*pi);
  Stat.ang      = prefAngle/2;
elseif strcmp(cfg.method,'direction')
  if (max(stimuli)-min(stimuli))>2*pi
    error('Directions has a range > 2*pi. Are you sure this is radians and not degrees?')
  end
  if (max(stimuli)-min(stimuli))<0.5*pi
    warning('Directions have a range < 1/2 pi. Are you sure this is correct?');
  end
  
  % compute the directionality index
  Stat.di = 1 - min(Tune.avg)./max(Tune.avg);
  
  % transform the data into complex numbers to compute resultant length
  z = exp(1i*stimuli(:)*ones(1,nUnits));
  sumAvg = sum(Tune.avg);
  z = Tune.avg.*z./(sumAvg(ones(nStimuli,1),:));
  Stat.osi  = abs(sum(z));
  
  % make preferred angle modulo 2pi, convert it back to range pi and convert to rad or deg
  prefAngle      = angle(nansum(z));
  Stat.ang = mod(prefAngle,2*pi);
end

Stat.label   = Tune.label;

% do the general cleanup and bookkeeping at the end of the function

ft_postamble previous   Tune
ft_postamble provenance Stat
ft_postamble history    Stat

