function [data] = denoise_synthetic(cfg, data);

% DENOISE_SYNTHETIC computes CTF higher-order synthetic gradients for
% preprocessed data and for the corresponding gradiometer definition.
%
% Use as
%   [data] = denoise_synthetic(cfg, data);
% where data should come from PREPROCESSING and the configuration should contain
%   cfg.gradient = 'none', 'G1BR', 'G2BR' or 'G3BR' specifies the gradiometer
%                  type to which the data should be changed
%   cfg.trials   = 'all' or a selection given as a 1xN vector (default = 'all')
%
% See also PREPROCESSING, DENOISE_SNS, DENOISE_TSR, DENOISE_PCA

% Copyright (C) 2004-2008, Robert Oostenveld
%
% $Log: denoise_synthetic.m,v $
% Revision 1.8  2009/07/08 07:20:48  roboos
% allow for arbitrary balancing and not only a small hardcoded list
%
% Revision 1.7  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.6  2008/06/26 16:09:15  roboos
% added trial selection
% added check on ctf input
% try to keep channel ordering identical
%
% Revision 1.5  2008/05/21 19:39:56  roboos
% small fix in documentation
%
% Revision 1.4  2008/05/21 10:24:41  roboos
% fixed input error for checkdata
%
% Revision 1.3  2008/05/20 11:54:07  roboos
% fixed typo, new should be data
%
% Revision 1.2  2008/05/15 15:07:46  roboos
% complete rewrite, initial implementation using apply_montage helper function
%
% Revision 1.1  2008/05/15 09:52:05  roboos
% renamed function syntheticgradient to denoise_synthetic
%

fieldtripdefs

data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

if ~senstype(data, 'ctf')
  error('synthetic gradients can only be computed for CTF data');
end

% set the defaults
if ~isfield(cfg, 'gradient'), error('cfg.gradient must be specified'); end
if ~isfield(cfg, 'trials'), cfg.trials = 'all'; end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
  data.offset = data.offset(cfg.trials);
  % update the trial definition (trl)
  if isfield(data, 'cfg') % try to locate the trl in the nested configuration
    trl = findcfg(data.cfg, 'trl');
  else
    trl = [];
  end
  if isempty(trl)
    % a trial definition is expected in each continuous data set
    warning('could not locate the trial definition ''trl'' in the data structure');
  else
    cfg.trlold=trl;
    cfg.trl=trl(cfg.trials,:);
  end
end

% remember the original channel ordering
labelorg = data.label;

% apply the balancing to the MEG data and to the gradiometer definition
current = data.grad.balance.current;
desired = cfg.gradient;

if ~strcmp(current, 'none')
  % first undo/invert the previously applied balancing
  try
    current_montage = getfield(data.grad.balance, data.grad.balance.current);
  catch
    error('unknown balancing for input data');
  end
  fprintf('converting from "%s" to "none"\n', current);
  data.grad = apply_montage(data.grad, current_montage, 'keepunused', 'yes', 'inverse', 'yes');
  data      = apply_montage(data     , current_montage, 'keepunused', 'yes', 'inverse', 'yes');
  data.grad.balance.current = 'none';
end % if

if ~strcmp(desired, 'none')
  % then apply the desired balancing
  try
    desired_montage = getfield(data.grad.balance, cfg.gradient);
  catch
    error('unknown balancing for input data');
  end
  fprintf('converting from "none" to "%s"\n', desired);
  data.grad = apply_montage(data.grad, desired_montage, 'keepunused', 'yes', 'inverse', 'yes');
  data      = apply_montage(data     , desired_montage, 'keepunused', 'yes', 'inverse', 'yes');
  data.grad.balance.current = desired;
end % if

% reorder the channels to stay close to the original ordering
[selorg, selnew] = match_str(labelorg, data.label);
if numel(selnew)==numel(labelorg)
  for i=1:numel(data.trial)
    data.trial{i} = data.trial{i}(selnew,:);
  end
  data.label = data.label(selnew);
else
  warning('channel ordering might have changed');
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: denoise_synthetic.m,v 1.8 2009/07/08 07:20:48 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
data.cfg = cfg;

