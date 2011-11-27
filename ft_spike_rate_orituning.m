function [Stat] = ft_spike_rate_orituning(cfg,Tune)

% FT_SPIKE_RATE_ORITUNING computes a model of the firing rate as a function
% of orientation or direction.
%
% Use as
%   [stat] = ft_spike_rate_tuning(cfg, tune)
%
% Input tune should be the output from FT_SPIKE_RATE_CONDITION
%
% Configurations:
%   cfg.stimuli  = should be an 1 x nConditions array of orientations or directions
%   cfg.method   = model to apply, implemented are 'orientation' and 'direction'

% Copyright (C) 2010, Martin Vinck
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'stimuli', 'method'});
cfg = ft_checkopt(cfg,'stimuli','doublevector');
cfg = ft_checkopt(cfg,'method', 'char', {'orientation', 'direction'});

% check whether trials were kept in the rate function
if ~isfield(Tune, 'avg'), error('MATLAB:ft_spike_rate_tuning:noFieldTrial',...
    'TUNE should contain the field avg');
end
stimuli = cfg.stimuli(:);
nConditions = size(Tune.avg,1);
if nConditions~=length(stimuli)
  error('MATLAB:ft_spike_rate_tuning:cfg:designWrongLength',...
    'Length of cfg.stimuli should match number of conditions in Tune.avg');
end

% get the unique stimuli and the number of directions, these should match
nStimuli    = length(stimuli);
nUnits      = size(Tune.avg,2);

% change the directions so it is a circular variable with range 2*pi
if strcmp(cfg.method,'orientation')
  if (max(stimuli)-min(stimuli))>pi
    error('MATLAB:ft_spike_rate_tuning:cfg.method',...
      'If cfg.tuningtype is "orientation", CFG.STIMULI should have range of pi');
  end
  stimuli = stimuli*2;
  
  if (max(stimuli)-min(stimuli))<0.5*pi
    warning('MATLAB:ft_spike_rate_tuning:orientationsSmallRange',...
      'Orientations have a range < 1/2 pi. Are you sure this is correct?');
  end
  
  % compute the directionality index
  Stat.directionIndex = 1 - min(Tune.avg)./max(Tune.avg);
  
  % transform the data into complex numbers to compute resultant length
  z = exp(1i*stimuli(:)*ones(1,nUnits));
  sumAvg = sum(Tune.avg);
  z = Tune.avg.*z./(sumAvg(ones(nStimuli,1),:));
  Stat.resLen  = abs(sum(z));
  
  % make preferred angle modulo 2pi, convert it back to range pi and convert to rad or deg
  prefAngle           = angle(nansum(z));
  prefAngle           = mod(prefAngle,2*pi);
  Stat.prefAngle      = prefAngle/2;
elseif strcmp(cfg.method,'direction')
  if (max(stimuli)-min(stimuli))>2*pi
    error('MATLAB:ft_spike_rate_tuning:directionsRangeTooLarge','%s\n%s',...
      'Directions has a range > 2*pi. Are you sure this is radians and not degrees?',...);
      'Please put the directions in a range of 2pi')
  end
  if (max(stimuli)-min(stimuli))<0.5*pi
    warning('MATLAB:ft_spike_rate_tuning:directionsSmallRange',...
      'Directions have a range < 1/2 pi. Are you sure this is correct?');
  end
  
  % compute the directionality index
  Stat.directionIndex = 1 - min(Tune.avg)./max(Tune.avg);
  
  % transform the data into complex numbers to compute resultant length
  z = exp(1i*stimuli(:)*ones(1,nUnits));
  sumAvg = sum(Tune.avg);
  z = Tune.avg.*z./(sumAvg(ones(nStimuli,1),:));
  Stat.resLen  = abs(sum(z));
  
  % make preferred angle modulo 2pi, convert it back to range pi and convert to rad or deg
  prefAngle      = angle(nansum(z));
  Stat.prefAngle = mod(prefAngle,2*pi);
end

Stat.label   = Tune.label;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous Tune
ft_postamble history Stat

