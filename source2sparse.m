function [source] = source2sparse(source);

% SOURCE2SPARCE removes the grid locations outside the brain from the source 
% reconstruction, thereby saving memory.
%
% This invalidates the fields that describe the grid, and also makes it
% more difficult to make a plot of each of the slices of the source volume.
% The original source structure can be recreated using SOURCE2FULL.
%
% Use as
%   [source] = source2sparse(source)
%
% See also SOURCE2FULL

% Copyright (C) 2004, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

if ~isfield(source, 'inside')
  warning('no gridpoints defined inside the brain');
  source.inside = [];
end

if ~isfield(source, 'outside')
  warning('no gridpoints defined outside the brain');
  source.outside = [];
end

inside  = source.inside;
outside = source.outside;

fprintf('total number of dipoles        : %d\n', length(inside)+length(outside));
fprintf('number of dipoles inside  brain: %d\n', length(inside));
fprintf('number of dipoles outside brain: %d\n', length(outside));

% first do the non-trial fields
[param]    = parameterselection('all', source);
trlparam   = strmatch('trial', param);
sel        = setdiff(1:length(param), trlparam);
param      = param(sel);

for j = 1:length(param)
  dat    = getsubfield(source, param{j});
  source = setsubfield(source, param{j}, dat(inside));
end

% then do the trial fields
if isfield(source, 'trial'),
  for j = 1:length(source.trial)
    tmpsource     = source.trial(j);
    tmpsource.dim = source.dim; % to fool parameterselection
    tmpparam      = parameterselection('all', tmpsource);
    for k = 1:length(tmpparam)
      dat       = getsubfield(tmpsource, tmpparam{k});
      tmpsource = setsubfield(tmpsource, tmpparam{k}, dat(inside));
    end
    tmpsource       = rmfield(tmpsource, 'dim');
    source.trial(j) = tmpsource;
  end
elseif isfield(source, 'trialA'),
  for j = 1:length(source.trialA)
    tmpsource     = source.trialA(j);
    tmpsource.dim = source.dim; % to fool parameterselection
    tmpparam      = parameterselection('all', tmpsource);
    for k = 1:length(tmpparam)
      dat       = getsubfield(tmpsource, tmpparam{k});
      tmpsource = setsubfield(tmpsource, tmpparam{k}, dat(inside));
    end
    tmpsource        = rmfield(tmpsource, 'dim');
    source.trialA(j) = tmpsource;
  end
elseif isfield(source, 'trialB'),
  for j = 1:length(source.trialB)
    tmpsource     = source.trialB(j);
    tmpsource.dim = source.dim; % to fool parameterselection
    tmpparam      = parameterselection('all', tmpsource);
    for k = 1:length(tmpparam)
      dat       = getsubfield(tmpsource, tmpparam{k});
      tmpsource = setsubfield(tmpsource, tmpparam{k}, dat(inside));
    end
    tmpsource        = rmfield(tmpsource, 'dim');
    source.trialB(j) = tmpsource;
  end
end

% update the inside, outside and source position
if isfield(source, 'inside')
  source.inside  = 1:length(inside);
end
if isfield(source, 'outside')
  source.outside = [];
end
if isfield(source, 'pos')
  source.pos     = source.pos(inside,:);
end
  
cfg = [];
% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = source.cfg; end
% remember the exact configuration details in the output 
source.cfg = cfg;

