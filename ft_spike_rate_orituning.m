function [Stat] = ft_spike_rate_orituning(cfg,Tune)

% FT_SPIKE_RATE_TUNING computes a model of the firing rate as a function of orientation or direction
% Input TUNE should be the output from FT_SPIKE_RATE_CONDITION
% Use as
%   [TUNE] = FT_SPIKE_RATE_TUNING(CFG,RATE)
%
% Configurations options (CFG):
%
%   cfg.stimuli  = should be an 1 x nConditions array of orientations or directions
%   cfg.method   = model to apply, implemented are 'orientation' and 'direction'. 
%

% Martin Vinck (C) 2010


% check whether trials were kept in the rate function

if ~isfield(Tune, 'avg'), error('MATLAB:ft_spike_rate_tuning:noFieldTrial',...
    'TUNE should contain the field avg'); 
end
if ~isfield(cfg,'stimuli'), error('MATLAB:ft_spike_rate_tuning:cfg:stimuliMissing'), end
if ~isfield(cfg,'method'), error('MATLAB:ft_spike_rate_tuning:cfg:stimuliMissing'),...
  'Please choose a method, "orientation", or "direction"', end
stimuli = cfg.stimuli(:);
if ~isrealvec(stimuli)
  error('MATLAB:ft_spike_rate_tuning:stimui',...
    'STIMULI should be a real vector');
end
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

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, k] = dbstack;
  cfg.version.name = st(k);
end
% remember the configuration details of the input data
try, cfg.previous = Tune.cfg; end
% remember the exact configuration details in the output 
Stat.cfg     = cfg;
