function [dipout] = music(dip, grad, headmodel, dat, varargin)

% MUSIC source localization using MUltiple SIgnal Classification
%
% This is a signal subspace method, which covers the techniques for
% multiple source localization by using the eigen structure of the
% measured data matrix.
%
% Use as
%   [dipout] = music(dip, grad, headmodel, dat, ...)
%
% Optional input arguments should be specified as key-value pairs and can be
%   'cov'              = data covariance matrix
%   'numcomponent'     = integer number
%   'feedback'         = 'none', 'gui', 'dial', 'textbar', 'text', 'textcr', 'textnl'
%   'reducerank'       = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%   'normalize'        = normalize the leadfield
%   'normalizeparam'   = parameter for depth normalization (default = 0.5)
%
% The original reference is
%   J.C. Mosher, P.S. Lewis and R.M. Leahy, "Multiple dipole modeling and
%   localization from spatiotemporal MEG data", IEEE Trans. Biomed.
%   Eng., pp 541-557, June, 1992.

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

% get the optional settings, or use the default value
cov            = ft_getopt(varargin, 'cov');
numcomponent   = ft_getopt(varargin, 'numcomponent');     % this is required, see below
feedback       = ft_getopt(varargin, 'feedback', 'text');
% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = ft_getopt(varargin, 'reducerank');
normalize      = ft_getopt(varargin, 'normalize');
normalizeparam = ft_getopt(varargin, 'normalizeparam');

if isempty(numcomponent)
  ft_error('you must specify the number of signal components');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside')
  dip.inside = ft_inside_vol(dip.pos, headmodel);
end

if any(dip.inside>1)
  % convert to logical representation
  tmp = false(size(dip.pos,1),1);
  tmp(dip.inside) = true;
  dip.inside = tmp;
end

% keep the original details on inside and outside positions
originside = dip.inside;
origpos    = dip.pos;

% select only the dipole positions inside the brain for scanning
dip.pos    = dip.pos(originside,:);
dip.inside = true(size(dip.pos,1),1);

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
dipout.jr = nan(size(dip.pos,1),1);

ft_progress('init', feedback, 'computing music metric');
for i=1:size(dip.pos,1)
  
  ft_progress(i/size(dip.pos,1), 'computing music metric %d/%d\n', i, size(dip.pos,1));
  
  if isfield(dip, 'leadfield')
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};
  elseif isfield(dip, 'mom')
    % compute the leadfield for a fixed dipole orientation
    lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(dip.pos(i,:), grad, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
  end
  
  % compute the MUSIC metric, c.f. equation 26
  dipout.jr(i) = (norm(ps * lf)./norm(lf)).^2;
  % as described in the Mosher 1992 paper on page 550, "...the general approach is to
  % evaluare Jr(i) over a fine three-dimensional grid, plot its inverse,
  % and look for p sharp spikes..."
  
end
ft_progress('close');

% wrap it all up, prepare the complete output
dipout.inside   = originside;
dipout.pos      = origpos;

% reassign the scan values over the inside and outside grid positions
if isfield(dipout, 'jr')
  dipout.jr( originside) = dipout.jr;
  dipout.jr(~originside) = nan;
end

