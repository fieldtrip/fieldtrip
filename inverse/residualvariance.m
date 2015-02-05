function [dipout] = residualvariance(dip, grad, vol, dat, varargin)

% RESIDUALVARIANCE scan with a single dipole and computes the RV
% at each grid location.
%
% Use as
%   [dipout] = residualvariance(dip, grad, vol, dat, ...)

% Copyright (C) 2004-2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% get the optional settings, or use default value
feedback      = keyval('feedback',      varargin); if isempty(feedback),      feedback = 'text';            end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside')
  dip.inside = ft_inside_vol(dip.pos, vol);
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
if isfield(dip, 'mom'),
  dip.mom = dip.mom(:,originside);
end
if isfield(dip, 'subspace')
  dip.subspace = dip.subspace(originside);
end

if isfield(dip, 'subspace')
  % remember the original data prior to the voxel dependant subspace projection
  dat_pre_subspace = dat;
  fprintf('using subspace projection\n');
end

mom = nan(size(dip.pos,1),1);
rv  = nan(size(dip.pos,1),1);
pow = nan(size(dip.pos,1),1);

ft_progress('init', feedback, 'computing inverse');
for i=1:size(dip.pos,1)
  
  ft_progress(i/size(dip.pos,1), 'computing inverse %d/%d\n', i, size(dip.pos,1));
  
  if isfield(dip, 'leadfield')
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};
  else
    % compute the leadfield
    lf = ft_compute_leadfield(dip.pos(i,:), grad, vol);
  end
  
  if isfield(dip, 'subspace')
    % do subspace projection of the forward model
    lf = dip.subspace{i} * lf;
    % the data and the covariance become voxel dependent due to the projection
    dat = dip.subspace{i} * dat_pre_subspace;
  end
  
  % compute spatiotemporal inverse using regional source
  lfi    = pinv(lf);
  mom{i} = lfi * dat;
  rv(i)  = sum(sum((dat - lf*mom{i}).^2, 1), 2)./sum(sum(dat.^2, 1), 2);
  
  % for plotting convenience also compute power at each location
  % FIXME is this normalization correct?
  pow(i) = mean(sum(mom{i}(:).^2, 1));
end
ft_progress('close');

% wrap it all up, prepare the complete output
dipout.inside  = originside;
dipout.pos     = origpos;

dipout.mom = nan(size(originside));
dipout.rv  = nan(size(originside));
dipout.pow = nan(size(originside));

% assign the output data
dipout.mom(originside) = mom(:);  % ensure that it is a column vector
dipout.rv(originside)  = rv(:);   % ensure that it is a column vector
dipout.pow(originside) = pow(:);  % ensure that it is a column vector
