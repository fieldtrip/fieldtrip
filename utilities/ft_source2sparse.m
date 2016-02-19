function [source] = ft_source2sparse(source)

% FT_SOURCE2SPARSE removes the grid locations outside the brain from the source 
% reconstruction, thereby saving memory.
%
% This invalidates the fields that describe the grid, and also makes it
% more difficult to make a plot of each of the slices of the source volume.
% The original source structure can be recreated using FT_SOURCE2FULL.
%
% Use as
%   [source] = ft_source2sparse(source)
%
% See also FT_SOURCE2FULL

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

ft_defaults

if ~isfield(source, 'inside')
  warning('no gridpoints defined inside the brain');
  source.inside = [];
elseif all(islogical(source.inside))
  source = fixinside(source, 'index'); % in contrast to the new convention, this function still relies on an indexed inside
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

% determine whether the source is old or new style
fnames = fieldnames(source);
if any(~cellfun('isempty', strfind(fnames, 'dimord'))),
  stype = 'new';
else
  stype = 'old';
end

if strcmp(stype, 'old'),
  % original code
  % first do the non-trial fields
  [param]    = parameterselection('all', source);
  %trlparam   = find(strcmp('trial', param));
  %sel        = setdiff(1:length(param), trlparam);
  %param      = param(sel);
  param      = setdiff(param, {'trial' 'pos'});
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
    source.inside  = [1:length(inside)]';
  end
  if isfield(source, 'outside')
    source.outside = [];
  end
  if isfield(source, 'pos')
    source.pos     = source.pos(inside,:);
  end
elseif strcmp(stype, 'new')
  % new style conversion
  nvox = numel(inside) + numel(outside);
  for k = 1:numel(fnames)
    tmpsiz = size(source.(fnames{k}));
    if any(tmpsiz==nvox)
      tmpsel = find(tmpsiz==nvox);
      if tmpsel==1,
        source.(fnames{k}) = source.(fnames{k})(inside,:,:,:,:);
      elseif tmpsel==2,
        source.(fnames{k}) = source.(fnames{k})(:,inside,:,:,:);
      else
        warning('not subselecting voxels, because location of pos-dimension is unexpected');
      end
    end
  end 
  
  % update the inside and outside
  if isfield(source, 'inside')
    source.inside  = [1:length(inside)]';
  end
  if isfield(source, 'outside')
    source.outside = [];
  end
end
  
cfg = [];
% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with MATLAB versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = source.cfg; end
% remember the exact configuration details in the output 
source.cfg = cfg;

