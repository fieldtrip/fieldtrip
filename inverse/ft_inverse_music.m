function [estimate] = ft_inverse_music(sourcemodel, sens, headmodel, dat, varargin)

% FT_INVERSE_MUSIC source localization using MUltiple SIgnal Classification.
% This is a signal subspace method, which covers the techniques for
% multiple source localization by using the eigen-structure of the
% measured data matrix.
%
% Use as
%   [estimate] = ft_inverse_music(sourcemodel, sens, headmodel, dat, ...)
% where
%   sourcemodel is the input source model, see FT_PREPARE_SOURCEMODEL
%   sens        is the gradiometer or electrode definition, see FT_DATATYPE_SENS
%   headmodel   is the volume conductor definition, see FT_PREPARE_HEADMODEL
%   dat         is the data matrix with the ERP or ERF
% and
%   estimate    contains the estimated source parameters
%
% Additional input arguments should be specified as key-value pairs and can include
%   'cov'              = data covariance matrix
%   'numcomponent'     = integer number
%   'feedback'         = can be 'none', 'gui', 'dial', 'textbar', 'text', 'textcr', 'textnl' (default = 'text')
%
% These options influence the forward computation of the leadfield
%   'reducerank'      = 'no' or number  (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no', in the case of a rank reduction this parameter determines whether the result will be backprojected onto the original subspace (default = 'yes')
%   'normalize'       = 'no', 'yes' or 'column' (default = 'no')
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate for the size of the corresponding patch (default = 1)
%
% This implements
% - J.C. Mosher, P.S. Lewis and R.M. Leahy, "Multiple dipole modeling and
%   localization from spatiotemporal MEG data", IEEE Trans. Biomed. Eng., 
%   pp 541-557, June, 1992.
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

% Copyright (C) 2004-2008, Robert Oostenveld
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

if mod(nargin-4,2)
  % the first 4 arguments are fixed, the other arguments should come in pairs
  ft_error('invalid number of optional arguments');
end

% get the optional input arguments, or use defaults
cov            = ft_getopt(varargin, 'cov');
numcomponent   = ft_getopt(varargin, 'numcomponent');     % this is required, see below
feedback       = ft_getopt(varargin, 'feedback', 'text');

% construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD
leadfieldopt = {};
leadfieldopt = ft_setopt(leadfieldopt, 'reducerank',     ft_getopt(varargin, 'reducerank'));
leadfieldopt = ft_setopt(leadfieldopt, 'backproject',    ft_getopt(varargin, 'backproject'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(varargin, 'normalize'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(varargin, 'normalizeparam'));
leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(varargin, 'weight'));

if isempty(numcomponent)
  ft_error('you must specify the number of signal components');
end

% flags to avoid calling isfield repeatedly in the loop over grid positions (saves a lot of time)
hasmom        = isfield(sourcemodel, 'mom');
hasleadfield  = isfield(sourcemodel, 'leadfield');
hasfilter     = false; % not used here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(sourcemodel, 'inside')
  if hasfilter
    sourcemodel.inside = ~cellfun(@isempty, sourcemodel.filter);
  elseif hasleadfield
    sourcemodel.inside = ~cellfun(@isempty, sourcemodel.leadfield);
  else
    sourcemodel.inside = ft_inside_headmodel(sourcemodel.pos, headmodel);
  end
end

% convert to logical representation
sourcemodel = fixinside(sourcemodel);

% keep the original details on inside and outside positions
originside = sourcemodel.inside;
origpos    = sourcemodel.pos;

% select only the dipole positions inside the brain for scanning
sourcemodel.pos    = sourcemodel.pos(originside,:);
sourcemodel.inside = true(size(sourcemodel.pos,1),1);

if hasmom
  sourcemodel.mom = sourcemodel.mom(:,originside);
end

if hasleadfield
  ft_info('using precomputed leadfields\n');
  sourcemodel.leadfield = sourcemodel.leadfield(originside);
else
  ft_info('computing forward model on the fly\n');
end

if ~isempty(cov)
  % compute signal and noise subspace from covariance matrix
  [u, s, v] = svd(cov);
else
  % compute signal and noise subspace from average data matrix
  [u, s, v] = svd(dat);
end
% select the noise subspace, c.f. equation 25
us = u(:,(numcomponent+1):end);
ps = us * us';

% allocate space to hold the result
estimate = [];
estimate.jr = nan(size(sourcemodel.pos,1),1);

ft_progress('init', feedback, 'scanning grid');
for i=1:size(sourcemodel.pos,1)
  ft_progress(i/size(sourcemodel.pos,1), 'scanning grid %d/%d\n', i, size(sourcemodel.pos,1));
  
  if hasleadfield && hasmom && size(sourcemodel.mom, 1)==size(sourcemodel.leadfield{i}, 2)
    % reuse the leadfield that was previously computed and project
    lf = sourcemodel.leadfield{i} * sourcemodel.mom(:,i);
  elseif  hasleadfield &&  hasmom
    % reuse the leadfield that was previously computed but don't project
    lf = sourcemodel.leadfield{i};
  elseif  hasleadfield && ~hasmom
    % reuse the leadfield that was previously computed
    lf = sourcemodel.leadfield{i};
  elseif ~hasleadfield &&  hasmom
    % compute the leadfield for a fixed dipole orientation
    lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:}) * sourcemodel.mom(:,i);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:});
  end
  
  % compute the music metric, c.f. equation 26
  estimate.jr(i) = (norm(ps * lf)./norm(lf)).^2;
  % as described in the Mosher 1992 paper on page 550, "...the general approach is to
  % evaluare Jr(i) over a fine three-dimensional grid, plot its inverse,
  % and look for p sharp spikes..."
  
end % for each dipole position
ft_progress('close');

% reassign the estimated values over the inside and outside grid positions
estimate.inside   = originside;
estimate.pos      = origpos;
if isfield(estimate, 'jr')
  estimate.jr( originside) = estimate.jr;
  estimate.jr(~originside) = nan;
end
if isfield(estimate, 'leadfield')
  estimate.leadfield( originside) = estimate.leadfield;
  estimate.leadfield(~originside) = {[]};
end
