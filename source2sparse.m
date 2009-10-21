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
% $Log: source2sparse.m,v $
% Revision 1.10  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.9  2006/03/30 12:24:33  roboos
% Implemented private/fixinside, which facilitates consistent
% handling of source/volume data. Improved documentation. Fixed some
% bugs related to inconsistent handling of ROIs (i.e. inside/outside)
%
% Revision 1.8  2006/01/31 12:57:20  jansch
% replaced explicit checking of all known parameters by parameterselection
%
% Revision 1.7  2005/08/23 13:15:16  jansch
% added source.avg.csd and noisecsd
%
% Revision 1.6  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.5  2005/01/25 11:08:23  roboos
% added var and sem as fields containing sourceparameters
%
% Revision 1.4  2004/09/23 15:07:16  roboos
% added support for lbex as cell-array in the source structure
%
% Revision 1.3  2004/08/26 12:15:25  roboos
% added some extra fprintf information
%
% Revision 1.2  2004/08/05 15:37:05  roboos
% added support for statistical parameters, fixed bug in leadfield
%
% Revision 1.1  2004/08/03 09:06:19  roboos
% initial implementation of these helper functions for beamformer sourceanalysis
%

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
cfg.version.id = '$Id: source2sparse.m,v 1.10 2008/09/22 20:17:44 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = source.cfg; end
% remember the exact configuration details in the output 
source.cfg = cfg;

