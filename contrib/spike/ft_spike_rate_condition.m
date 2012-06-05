function [Tune] = ft_spike_rate_condition(cfg,Rate)

% FT_SPIKE_RATE_CONDITION computes the average rate for different
% conditions. This function should be the precursor for an orientation
% tuning / constrast tuning script.
%
% Use as
%   [tune] = ft_spike_rate_condition(cfg, rate)
%
% The input variable rate should be the output from FT_SPIKE_RATE
% with cfg.keeptrials = 'yes'.
%
% Configurations:
%   cfg.design = should be an 1 x nTrials array, with an integer
%                value for every condition
%
% See also FT_SPIKE_RATE

% Martin Vinck (C) 2010
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'design'});
cfg = ft_checkopt(cfg,'design','doublevector');

% check whether trials were kept in the rate function
if ~isfield(Rate, 'trial'), error('MATLAB:ft_spike_rate_condition:noFieldTrial',...
    'RATE should contain the field trial (use cfg.keeptrials = "yes" in spike_rate)');
end
design = cfg.design(:);

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

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous Rate
ft_postamble history Tune

