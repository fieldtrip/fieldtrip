function [Tune] = ft_spike_rate_condition(cfg,Rate)

% FT_SPIKE_RATE_CONDITION computes the average rate for different conditions.
% Input RATE should be the output from FT_SPIKE_RATE with cfg.keeptrials = 'yes'
% Use as
%   [TUNE] = FT_SPIKE_RATE_CONDITION(CFG,RATE)
%
% Configurations options (CFG):
%
%   cfg.design      = should be an 1 x nTrials array, with an integer value for every condition
%
% This function should be the precursor for an orientation tuning / constrast tuning script

% Martin Vinck (C) 2010

% check whether trials were kept in the rate function
if ~isfield(Rate, 'trial'), error('MATLAB:ft_spike_rate_condition:noFieldTrial',...
    'RATE should contain the field trial (use cfg.keeptrials = "yes" in spike_rate)'); 
end
if ~isfield(cfg,'design'), error('MATLAB:ft_spike_rate_condition:cfg:designMissing'), end
design = cfg.design(:);
if ~isrealvec(design)
  error('MATLAB:ft_spike_rate_condition:design',...
    'DESIGN should be a real vector');
end
nTrials = size(Rate.trial,1);
if nTrials~=length(design)
  error('MATLAB:ft_spike_rate_condition:cfg:designWrongLength',...
    'Length of cfg.design should match number of trials in RATE.trial');
end

% get the unique stimuli and the number of directions, these should match
conditions     = unique(design);
nConditions    = length(conditions);
if  nConditions==1
     error('MATLAB:ft_spike_rate_condition:numUniqueStim',...
    'number of unique elements in cfg.design should be >1 for this function to have meaning')
end
    
% compute the firing rate per stimulus condition
nUnits = length(Rate.avg);
[avg,var] = deal(NaN(nConditions,nUnits));
dof = zeros(nConditions,1);

for iCondition = 1:nConditions   
    hasStim            = design==conditions(iCondition);       
    dof(iCondition)    = sum(hasStim);
    avg(iCondition,:)  = nanmean(Rate.trial(hasStim,:),1); % 1 was missing! bug with 1 unit
    var(iCondition,:)  = nanvar(Rate.trial(hasStim,:),[],1);
end

% collect the results
Tune.avg            = avg;
Tune.var            = var;
Tune.dof            = dof;
Tune.label          = Rate.label;

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
try, cfg.previous = Rate.cfg; end
% remember the exact configuration details in the output 
Tune.cfg     = cfg;
